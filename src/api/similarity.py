# api/similarity.py

from fastapi import APIRouter, Query
from rdkit import Chem
import pandas as pd

from src.services.chemdb_service import get_engine
from data_sources import chembl_client
from src.services.conformer_manager import get_or_build_conformer
from src.services.pipeline import (
    compute_2d_similarity,
    compute_usr_similarity,
)

router = APIRouter(prefix="/similarity", tags=["similarity"])


@router.get("/")
def similarity_with_local_drugs(
    drug: str = Query(...),
    alpha: float = 0.7,
    radius: int = 2,
    nbits: int = 2048,
):
    """
    Computes similarity between the query drug and ALL local drugs stored in local_drugs.
    """

    # Step 1: resolve query drug via ChEMBL
    ref_df = chembl_client.search_molecule(drug)
    if len(ref_df) == 0:
        return {"error": f"{drug} not found in ChEMBL"}

    row = ref_df.iloc[0]
    ref_smiles = row["smiles"]
    ref_chembl_id = row["chembl_id"]

    ref_mol, source = get_or_build_conformer(
        normalized_name=f"_query_{drug.lower()}",
        chembl_id=ref_chembl_id,
        chembl_name=drug,
        smiles=ref_smiles,
        store_in_db=False,
    )

    if ref_mol is None:
        return {"error": "Could not build/load conformer for query drug"}

    ref_conf_ids = [conf.GetId() for conf in ref_mol.GetConformers()]

    # Step 2: Load local drugs
    engine = get_engine()
    df = pd.read_sql(
        "SELECT normalized_name, chembl_name, smiles, molblock FROM local_drugs",
        engine,
    )

    names = df["normalized_name"].tolist()
    smiles_list = df["smiles"].tolist()
    molblocks = df["molblock"].tolist()

    # Step 3: Build RDKit mols
    target_mols = []
    target_conf_ids = []

    for block in molblocks:
        mol = Chem.MolFromMolBlock(block, removeHs=False)
        if mol:
            target_mols.append(mol)
            target_conf_ids.append([conf.GetId() for conf in mol.GetConformers()])
        else:
            target_mols.append(None)
            target_conf_ids.append([])

    # Step 4: compute similarities
    df2d = compute_2d_similarity(ref_smiles, smiles_list, names, radius, nbits)
    df3d = compute_usr_similarity(ref_mol, ref_conf_ids, target_mols, target_conf_ids, names)

    merged = df2d.merge(df3d, on="name")
    merged["weighted_similarity"] = (
        alpha * merged["3D_similarity"]
        + (1 - alpha) * merged["2D_similarity"]
    )
    merged = merged.sort_values("weighted_similarity", ascending=False)

    return {
        "query": drug,
        "count": len(merged),
        "results": merged.to_dict(orient="records"),
    }
