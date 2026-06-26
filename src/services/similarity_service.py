import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdMolDescriptors
from sqlalchemy import text

from src.services.chemdb_service import get_engine
from src.data_sources import chembl_client
from src.services.conformer_manager import build_query_conformer


def _get_selected_target_ids(engine):
    query = text("""
        SELECT target_chembl_id
        FROM similarity_selected_targets
        WHERE is_active = TRUE
    """)

    with engine.connect() as conn:
        rows = conn.execute(query).fetchall()

    return {row[0] for row in rows}


def _get_drug_targets_map(engine, chembl_ids: list[str]):
    if not chembl_ids:
        return {}

    query = text("""
        SELECT
            molecule_chembl_id,
            target_chembl_id,
            target_name
        FROM drug_targets
        WHERE molecule_chembl_id IN :chembl_ids
    """)

    with engine.connect() as conn:
        rows = conn.execute(
            query.bindparams(chembl_ids=tuple(chembl_ids))
        ).fetchall()

    result = {}

    for molecule_chembl_id, target_chembl_id, target_name in rows:
        result.setdefault(molecule_chembl_id, {})
        result[molecule_chembl_id][target_chembl_id] = target_name or target_chembl_id

    return result


def enrich_with_shared_selected_targets(
    merged: pd.DataFrame,
    ref_chembl_id: str,
    engine,
):
    selected_target_ids = _get_selected_target_ids(engine)

    if not selected_target_ids:
        merged["shared_selected_targets"] = ""
        merged["shared_selected_target_count"] = 0
        return merged

    result_chembl_ids = merged["chembl_id"].dropna().tolist()
    all_chembl_ids = [ref_chembl_id] + result_chembl_ids

    targets_map = _get_drug_targets_map(engine, all_chembl_ids)

    ref_targets = targets_map.get(ref_chembl_id, {})
    ref_selected_target_ids = set(ref_targets.keys()) & selected_target_ids

    shared_names = []
    shared_counts = []

    for chembl_id in merged["chembl_id"]:
        candidate_targets = targets_map.get(chembl_id, {})
        shared_ids = ref_selected_target_ids & set(candidate_targets.keys())

        names = [
            ref_targets.get(target_id)
            or candidate_targets.get(target_id)
            or target_id
            for target_id in sorted(shared_ids)
        ]

        shared_names.append(", ".join(names))
        shared_counts.append(len(shared_ids))

    merged["shared_selected_targets"] = shared_names
    merged["shared_selected_target_count"] = shared_counts

    return merged

# ----------------------------
# 3D Similarity (multi confs)
# ----------------------------
def compute_usr_similarity(ref_mol, ref_confs, t_mols, t_confs, names):
    ref_fps = [
        np.array(rdMolDescriptors.GetUSRCAT(ref_mol, confId=cid))
        for cid in ref_confs
    ]

    results = []

    for mol, conf_ids, name in zip(t_mols, t_confs, names):
        if mol is None or not conf_ids:
            results.append((name, 0.0))
            continue

        best = 0.0
        for ref_fp in ref_fps:
            ref_norm = np.linalg.norm(ref_fp)
            if ref_norm == 0:
                continue

            for cid in conf_ids:
                try:
                    fp = np.array(rdMolDescriptors.GetUSRCAT(mol, confId=cid))
                    norm = np.linalg.norm(fp)
                    if norm == 0:
                        continue
                    sim = float(np.dot(ref_fp, fp) / (ref_norm * norm))
                    best = max(best, sim)
                except Exception:
                    continue

        results.append((name, best))

    return pd.DataFrame(results, columns=["name", "3D_similarity"])


# ----------------------------
# 2D Similarity
# ----------------------------
def compute_2d_similarity(ref_smiles, smiles_list, names, radius=2, nbits=2048):
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, radius, nBits=nbits)

    rows = []
    for s, name in zip(smiles_list, names):
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            rows.append((name, 0.0))
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
        sim = DataStructs.TanimotoSimilarity(ref_fp, fp)
        rows.append((name, sim))

    return pd.DataFrame(rows, columns=["name", "2D_similarity"])


def compute_similarity_with_local_drugs(
    drug: str,
    alpha: float,
    radius: int,
    nbits: int,
):
    # ------------------------------------------------------------
    # 1) Resolve query drug via ChEMBL
    # ------------------------------------------------------------
    query_norm = drug.strip().lower()
    ref_df = chembl_client.search_molecule(drug)
    if ref_df.empty:
        raise ValueError(f"Query drug '{drug}' not found in ChEMBL")

    ref_row = ref_df.iloc[0]
    ref_smiles = ref_row["smiles"]

    # ------------------------------------------------------------
    # 2) Build query conformer (temporary)
    # ------------------------------------------------------------
    ref_mol = build_query_conformer(ref_smiles)
    if ref_mol is None or ref_mol.GetNumConformers() == 0:
        raise RuntimeError("Failed to build 3D conformer for query drug")

    ref_conf_ids = [conf.GetId() for conf in ref_mol.GetConformers()]

    # ------------------------------------------------------------
    # 3) Load local drugs from DB
    # ------------------------------------------------------------
    engine = get_engine()
    df = pd.read_sql(
        """
        SELECT
            molecule_chembl_id,
            normalized_name,
            smiles,
            molblock
        FROM local_drugs
        WHERE molblock IS NOT NULL
        """,
        engine,
    )
    df = df[df["normalized_name"] != query_norm]

    if df.empty:
        raise RuntimeError("Local drug database is empty")
    
    chembl_ids = df["molecule_chembl_id"].tolist()
    names = df["normalized_name"].tolist()
    smiles_list = df["smiles"].tolist()
    molblocks = df["molblock"].tolist()

    # ------------------------------------------------------------
    # 4) Build RDKit mols from MolBlocks
    # ------------------------------------------------------------
    target_mols = []
    target_conf_ids = []

    for block in molblocks:
        mol = Chem.MolFromMolBlock(block, removeHs=False)
        if mol and mol.GetNumConformers() > 0:
            target_mols.append(mol)
            target_conf_ids.append([conf.GetId() for conf in mol.GetConformers()])
        else:
            target_mols.append(None)
            target_conf_ids.append([])

    # ------------------------------------------------------------
    # 5) Compute similarities
    # ------------------------------------------------------------
    df2d = compute_2d_similarity(
        ref_smiles=ref_smiles,
        smiles_list=smiles_list,
        names=names,
        radius=radius,
        nbits=nbits,
    )

    df3d = compute_usr_similarity(
        ref_mol=ref_mol,
        ref_confs=ref_conf_ids,
        t_mols=target_mols,
        t_confs=target_conf_ids,
        names=names,
    )

    merged = df2d.merge(df3d, on="name", how="inner")
    name_to_chembl = dict(zip(names, chembl_ids))
    merged["chembl_id"] = merged["name"].map(name_to_chembl)


    ref_chembl_id = ref_row.get("chembl_id") or ref_row.get("molecule_chembl_id")

    if ref_chembl_id:
        merged = enrich_with_shared_selected_targets(
            merged=merged,
            ref_chembl_id=ref_chembl_id,
            engine=engine,
        )
    else:
        merged["shared_selected_targets"] = ""
        merged["shared_selected_target_count"] = 0



    merged["weighted_similarity"] = (
        alpha * merged["3D_similarity"]
        + (1.0 - alpha) * merged["2D_similarity"]
    )

    merged = merged.sort_values(
        ["weighted_similarity", "shared_selected_target_count"],
        ascending=[False, False],
    ).reset_index(drop=True)

    cols_to_round = ["2D_similarity", "3D_similarity", "weighted_similarity"]
    merged[cols_to_round] = merged[cols_to_round].round(3)

    return merged

