# src/api/import_local_drugs.py

from fastapi import APIRouter, UploadFile, File, Query

from src.data_sources.chembl_client import get_local_approved_drugs
from src.services.drug_name_normalizer import normalize_raw_names
from src.services.chembl_lookup_service import build_chembl_lookup
from src.services.local_drug_service import (
    load_local_drugs_from_excel,
    process_local_drugs,
)

router = APIRouter(prefix="/import", tags=["import-local-drugs"])


@router.post("/local_drugs")
async def import_local_drugs(
    file: UploadFile = File(...),
    chembl_limit: int = Query(5000),
):
    print("\n================ Importing Local Drug List ================")
    content = await file.read()
    df = load_local_drugs_from_excel(content)

    # STEP 1: normalize local names
    all_norms = normalize_raw_names(df["Name"].astype(str).tolist())
    print(f"[STEP-1] Normalized names: {len(all_norms)}")

    # STEP 2: load approved drugs from local ChEMBL cache
    drugs_df = get_local_approved_drugs(limit=chembl_limit)
    print(
        f"[STEP-2] Loaded {len(drugs_df)} approved drugs from local ChEMBL cache "
        f"(limit={chembl_limit})"
    )

    # STEP 3: build normalized lookup for ChEMBL names
    chembl_lookup = build_chembl_lookup(drugs_df)
    print(f"[STEP-3] Normalized ChEMBL names in cache: {len(chembl_lookup)}")

    # STEP 4: process local drugs (DB + conformers)
    print("[STEP-4] Processing normalized local drug names...")
    inserted_count, unmatched = process_local_drugs(all_norms, chembl_lookup)

    print(f"[STEP-4] Inserted new drugs: {inserted_count}")
    print("================ Import Finished ================\n")

    return {
        "processed": len(all_norms),
        "inserted": inserted_count,
        "unmatched": unmatched,
    }
