# src/services/local_drug_service.py

import io
from typing import Dict, List, Tuple

import pandas as pd
from rdkit import Chem

from src.services.conformer_manager import get_or_build_conformer
from src.services.chemdb_service import get_local_drug_by_normalized_name


def load_local_drugs_from_excel(content: bytes) -> pd.DataFrame:
    """
    Load Excel drug list.
    A column named 'Name' is expected.
    """
    buf = io.BytesIO(content)
    df = pd.read_excel(buf)

    if "Name" not in df.columns:
        raise ValueError("Column 'Name' not found in Excel file.")

    df = df[df["Name"].notna()].copy()
    return df


def process_local_drugs(
    all_norms: List[str],
    chembl_lookup: Dict[str, pd.Series],
) -> Tuple[int, List[str]]:
    """
    Process normalized local drug names:
      - Skip local entries that already exist
      - Match against local ChEMBL cache
      - Validate SMILES
      - Reject large molecules
      - Build or load conformer

    Args:
        all_norms: normalized local drug names
        chembl_lookup: mapping normalized_name -> ChEMBL row (pandas Series)

    Returns:
        inserted_count: number of successfully inserted drugs
        unmatched_list: list of names that could not be processed
    """

    inserted_count = 0
    unmatched: List[str] = []

    for norm in all_norms:
        existing = get_local_drug_by_normalized_name(norm)
        if existing:
            print(f"  ✓ {norm}: already exists")
            continue

        if norm not in chembl_lookup:
            print(f"  ✗ {norm}: not found in local ChEMBL cache → skip")
            unmatched.append(norm)
            continue

        row = chembl_lookup[norm]
        chembl_id = row["molecule_chembl_id"]
        smiles = row["smiles"]
        chembl_name = row["name"]

        print(f"  → New drug: {norm} | chembl={chembl_id}")

        # Validate SMILES
        mol0 = Chem.MolFromSmiles(smiles)
        if mol0 is None:
            print("     ✗ Invalid SMILES → skip")
            unmatched.append(norm)
            continue

        # Reject molecules that are too large
        if mol0.GetNumAtoms() > 120:
            print("     ✗ Molecule too large → skip")
            unmatched.append(norm)
            continue

        # Build or load conformer
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
            unmatched.append(norm)

    return inserted_count, unmatched
