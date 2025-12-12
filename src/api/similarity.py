# src/api/similarity.py

from fastapi import APIRouter, Query
from rdkit import Chem
import pandas as pd

from src.services.chemdb_service import get_engine
from src.data_sources import chembl_client
from src.services.conformer_manager import build_query_conformer
from src.services.pipeline import (
    compute_2d_similarity,
    compute_usr_similarity,
)

router = APIRouter(prefix="/similarity", tags=["similarity"])


@router.get("/")
def similarity_with_local_drugs(
    drug: str = Query(..., description="Query drug name"),
    alpha: float = Query(0.7, ge=0.0, le=1.0),
    radius: int = Query(2),
    nbits: int = Query(2048),
):
    """
    Compute similarity between a query drug and ALL locally stored drugs.

    Similarity = alpha * 3D_USRCAT + (1 - alpha) * 2D_Tanimoto
    """

    # ------------------------------------------------------------
    # STEP 1: Resolve query drug via ChEMBL (SMILES only)
    # ------------------------------------------------------------
    ref_df = chembl_client.search_molecule(drug)
    if ref_df.empty:
        return {"error": f"Query drug '{drug}' not found in ChEMBL"}

    ref_row = ref_df.iloc[0]
    ref_smiles = ref_row["smiles"]
    ref_chembl_id = ref_row["chembl_id"]

    # ------------------------------------------------------------
    # STEP 2: Build query conformer (NOT stored in DB)
    # ------------------------------------------------------------
    ref_mol = build_query_conformer(
        chembl_id=ref_chembl_id,
        smiles=ref_smiles,
    )

    if ref_mol is None or ref_mol.GetNumConformers() == 0:
        return {"error": "Failed to build 3D conformer for query drug"}

    ref_conf_ids = [conf.GetId() for conf in ref_mol.GetConformers()]

    # ------------------------------------------------------------
    # STEP 3: Load local drug dataset from DB
    # ------------------------------------------------------------
    engine = get_engine()
    df = pd.read_sql(
        """
        SELECT
            normalized_name,
            chembl_name,
            smiles,
            molblock
        FROM local_drugs
        WHERE molblock IS NOT NULL
        """,
        engine,
    )

    if df.empty:
        return {"error": "Local drug database is empty"}

    names = df["normalized_name"].tolist()
    smiles_list = df["smiles"].tolist()
    molblocks = df["molblock"].tolist()

    # ------------------------------------------------------------
    # STEP 4: Build RDKit molecules from stored MolBlocks
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
    # STEP 5: Compute similarities
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
        ref_conf_ids=ref_conf_ids,
        target_mols=target_mols,
        target_conf_ids=target_conf_ids,
        names=names,
    )

    merged = df2d.merge(df3d, on="name", how="inner")

    merged["weighted_similarity"] = (
        alpha * merged["3D_similarity"]
        + (1.0 - alpha) * merged["2D_similarity"]
    )

    merged = merged.sort_values(
        "weighted_similarity", ascending=False
    ).reset_index(drop=True)

    # ------------------------------------------------------------
    # STEP 6: Return response
    # ------------------------------------------------------------
    return {
        "query": drug,
        "alpha": alpha,
        "count": len(merged),
        "results": merged.to_dict(orient="records"),
    }
