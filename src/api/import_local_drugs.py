# api/import_local_drugs.py

from fastapi import APIRouter, UploadFile, File
from rdkit import Chem

from data_sources.chembl_client import get_local_approved_drugs
from src.services.name_normalizer import normalize_raw_names
from src.services.chembl_lookup_service import build_chembl_lookup
from src.services.chemdb_service import get_local_drug_by_normalized_name
from src.services.conformer_manager import get_or_build_conformer
from src.services.local_drug_service import load_local_drugs_from_excel

router = APIRouter(prefix="/import", tags=["import-local-drugs"])


@router.post("/local_drugs")
async def import_local_drugs(
    file: UploadFile = File(...),
):
    print("\n================ Importing Local Drug List ================")
    content = await file.read()
    df = load_local_drugs_from_excel(content)

    # STEP 1 — Normalize names
    all_norms = normalize_raw_names(df["Name"].astype(str).tolist())
    print(f"[STEP-1] Normalized names: {len(all_norms)}")

    # STEP 2 — Load approved ChEMBL drugs from LOCAL database
    drugs_df = get_local_approved_drugs()
    print(f"[STEP-2] Loaded {len(drugs_df)} approved drugs from local DB")

    # STEP 3 — Build lookup from local chembl_approved_drugs
    chembl_lookup = build_chembl_lookup(drugs_df)
    print(f"[STEP-3] Normalized ChEMBL names (local): {len(chembl_lookup)}")

    print("[STEP-4] Processing normalized local drug names...")

    inserted_count = 0

    # STEP 4 — Incremental processing
    for norm in all_norms:
        existing = get_local_drug_by_normalized_name(norm)
        if existing:
            print(f"  ✓ {norm}: already exists")
            continue

        if norm not in chembl_lookup:
            print(f"  ✗ {norm}: not found in local ChEMBL cache → skip")
            continue

        row = chembl_lookup[norm]
        chembl_id = row["molecule_chembl_id"]
        smiles = row["smiles"]
        chembl_name = row["name"]

        print(f"  → New drug: {norm} | chembl={chembl_id}")

        mol0 = Chem.MolFromSmiles(smiles)
        if mol0 is None:
            print("     ✗ Invalid SMILES → skip")
            continue

        if mol0.GetNumAtoms() > 120:
            print("     ✗ Molecule too large → skip")
            continue

        # generate or load conformer
        mol, source = get_or_build_conformer(
            normalized_name=norm,
            chembl_id=chembl_id,
            chembl_name=chembl_name,
            smiles=smiles,
        )

        if mol:
            print(f"     ✓ Conformer ready | from: {source}")
            inserted_count += 1
        else:
            print("     ✗ Conformer generation failed → skip")

    print(f"[STEP-4] Inserted new drugs: {inserted_count}")
    print("================ Import Finished ================\n")

    return {
        "processed": len(all_norms),
        "inserted": inserted_count,
    }

