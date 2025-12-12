# conformer_manager.py

from typing import Optional, Tuple
from rdkit import Chem
from src.services.chemdb_service import (
    get_local_drug_by_normalized_name,
    bulk_upsert_local_drugs,
)
from src.services.pubchem_utils import load_pubchem_conformer_from_smiles


def mol_from_molblock(block: str) -> Optional[Chem.Mol]:
    """
    Convert MolBlock text to an RDKit Mol object.
    """
    if not block:
        return None
    return Chem.MolFromMolBlock(block, removeHs=False)


def mol_to_molblock(mol: Chem.Mol) -> str:
    """
    Convert RDKit Mol to V2000 MolBlock text.
    """
    return Chem.MolToMolBlock(mol, forceV3000=False)


def generate_rdkit_conformer(smiles: str) -> Optional[Chem.Mol]:
    """
    Build a 3D conformer using RDKit ETKDG.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)

    from rdkit.Chem import AllChem
    params = AllChem.ETKDGv3()
    params.randomSeed = 42

    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=1, params=params)
    if len(conf_ids) == 0:
        return None

    AllChem.UFFOptimizeMolecule(mol)
    return mol


def get_or_build_conformer(
    normalized_name: str,
    chembl_id: str,
    chembl_name: str,
    smiles: str,
):
    """
    Load stored conformer if exists, otherwise build it from:
      1) PubChem bioactive conformer
      2) RDKit fallback
    Store the result inside local_drugs table.
    """

    # 1) Try loading from DB
    row = get_local_drug_by_normalized_name(normalized_name)
    if row and row.get("molblock"):
        mol = mol_from_molblock(row["molblock"])
        if mol:
            return mol, row.get("conformer_source")

    # 2) Try PubChem
    mol = load_pubchem_conformer_from_smiles(smiles)
    if mol:
        molblock = mol_to_molblock(mol)

        bulk_upsert_local_drugs([{
            "normalized_name": normalized_name,
            "local_name": normalized_name,
            "chembl_name": chembl_name,
            "molecule_chembl_id": chembl_id,
            "smiles": smiles,
            "conformer_source": "pubchem",
            "molblock": molblock,
        }])

        return mol, "pubchem"

    # 3) RDKit fallback
    mol = generate_rdkit_conformer(smiles)
    if mol:
        molblock = mol_to_molblock(mol)

        bulk_upsert_local_drugs([{
            "normalized_name": normalized_name,
            "local_name": normalized_name,
            "chembl_name": chembl_name,
            "molecule_chembl_id": chembl_id,
            "smiles": smiles,
            "conformer_source": "rdkit",
            "molblock": molblock,
        }])

        return mol, "rdkit"

    return None, None
