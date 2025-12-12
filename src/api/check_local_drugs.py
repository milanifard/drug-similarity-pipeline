# src/api/check_local_drugs.py

from fastapi import APIRouter, UploadFile, File
from services.name_normalizer import normalize_raw_names
from services.local_drug_service import load_local_drugs_from_excel
from services.chembl_lookup_service import build_chembl_lookup
from src.data_sources.chembl_client import get_local_approved_drugs

router = APIRouter(prefix="/check", tags=["drug-check"])

@router.post("/local_drugs")
async def check_local_drugs(
    file: UploadFile = File(...)
):
    """
    Reads the uploaded Excel file, normalizes drug names,
    and checks which ones can be matched against the local ChEMBL cache.

    It does NOT perform conformer generation.
    """

    # Load excel
    content = await file.read()
    df = load_local_drugs_from_excel(content)

    # Normalize names
    all_norms = normalize_raw_names(df["Name"].astype(str).tolist())

    # Load local ChEMBL cache (already synced)
    chembl_df = get_local_approved_drugs(limit=50000)
    chembl_lookup = build_chembl_lookup(chembl_df)

    matched = []
    unmatched = []

    for norm in all_norms:
        if norm in chembl_lookup:
            matched.append(norm)
        else:
            unmatched.append(norm)

    return {
        "total_normalized": len(all_norms),
        "matched": matched,
        "matched_count": len(matched),
        "unmatched": unmatched,
        "unmatched_count": len(unmatched),
    }
