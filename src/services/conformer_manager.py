# src/services/conformer_manager.py

from typing import Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem

from src.services.chemdb_service import (
    get_local_drug_by_normalized_name,
    bulk_upsert_local_drugs,
)
from src.services.pubchem_utils import load_pubchem_conformer_from_smiles


# ============================================================
# Core utilities (NO DB, NO POLICY)
# ============================================================

def mol_from_molblock(block: str) -> Optional[Chem.Mol]:
    """Convert MolBlock text to RDKit Mol."""
    if not block:
        return None
    return Chem.MolFromMolBlock(block, removeHs=False)


def mol_to_molblock(mol: Chem.Mol) -> str:
    """Convert RDKit Mol to V2000 MolBlock."""
    return Chem.MolToMolBlock(mol, forceV3000=False)


def build_conformer_from_smiles(
    smiles: str,
    num_confs: int,
) -> Optional[Chem.Mol]:
    """
    Build 3D conformers using RDKit ETKDG.

    This function:
      - does NOT touch DB
      - does NOT decide policy
      - only builds conformers
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.maxAttempts = 5000
    params.pruneRmsThresh = 0.1

    conf_ids = AllChem.EmbedMultipleConfs(
        mol,
        numConfs=num_confs,
        params=params,
    )

    if not conf_ids:
        return None

    # Optimize all conformers
    AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=600)
    return mol


# ============================================================
# PubChem strategy (isolated)
# ============================================================

def try_pubchem_conformer(smiles: str) -> Optional[Chem.Mol]:
    """
    Try to load a bioactive 3D conformer from PubChem.
    Returns Mol if available, otherwise None.
    """
    mol = load_pubchem_conformer_from_smiles(smiles)
    if mol and mol.GetNumConformers() > 0:
        AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=600)
        return mol
    return None


# ============================================================
# Query conformer (NO persistence)
# ============================================================

def build_query_conformer(smiles: str) -> Optional[Chem.Mol]:
    """
    Build 3D conformer for QUERY molecule.

    Priority:
      1) PubChem bioactive conformer
      2) RDKit ETKDG (multiple conformers)

    This function NEVER writes to DB.
    """

    # 1) PubChem
    mol = try_pubchem_conformer(smiles)
    if mol:
        return mol

    # 2) RDKit fallback (higher accuracy)
    return build_conformer_from_smiles(
        smiles=smiles,
        num_confs=30,
    )


# ============================================================
# Local drug conformer (WITH persistence)
# ============================================================

def get_or_build_local_conformer(
    normalized_name: str,
    chembl_id: str,
    chembl_name: str,
    smiles: str,
) -> Tuple[Optional[Chem.Mol], Optional[str]]:
    """
    Load or build conformer for LOCAL drug.

    Priority:
      1) Load from local DB
      2) PubChem bioactive conformer
      3) RDKit ETKDG fallback

    The result IS stored in local_drugs table.
    """

    # --------------------------------------------------------
    # 1) Load from DB
    # --------------------------------------------------------
    row = get_local_drug_by_normalized_name(normalized_name)
    if row and row.get("molblock"):
        mol = mol_from_molblock(row["molblock"])
        if mol:
            return mol, row.get("conformer_source")

    # --------------------------------------------------------
    # 2) Try PubChem
    # --------------------------------------------------------
    mol = try_pubchem_conformer(smiles)
    source = "pubchem"

    # --------------------------------------------------------
    # 3) RDKit fallback
    # --------------------------------------------------------
    if not mol:
        mol = build_conformer_from_smiles(
            smiles=smiles,
            num_confs=30,
        )
        source = "rdkit"

    # --------------------------------------------------------
    # 4) Persist if successful
    # --------------------------------------------------------
    if mol:
        bulk_upsert_local_drugs([{
            "normalized_name": normalized_name,
            "local_name": normalized_name,
            "chembl_name": chembl_name,
            "molecule_chembl_id": chembl_id,
            "smiles": smiles,
            "conformer_source": source,
            "molblock": mol_to_molblock(mol),
        }])
        return mol, source

    return None, None
