import requests
from rdkit import Chem

def get_pubchem_cid_from_smiles(smiles):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
    try:
        data = requests.get(url, timeout=10).json()
        return data["IdentifierList"]["CID"][0]
    except Exception:
        return None

def load_pubchem_conformer_from_smiles(smiles):
    cid = get_pubchem_cid_from_smiles(smiles)
    if cid is None:
        return None
    return load_pubchem_conformer_by_cid(cid)


def get_pubchem_cid(name):
    """Return PubChem CID for a given molecule name."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
    try:
        data = requests.get(url, timeout=10).json()
        return data["IdentifierList"]["CID"][0]
    except Exception:
        return None


def load_pubchem_conformer_by_cid(cid):
    """Download BioActive 3D conformer SDF from PubChem"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{cid}/record/SDF/?record_type=3d"
    try:
        sdf = requests.get(url, timeout=10).text
        mol = Chem.MolFromMolBlock(sdf, removeHs=False)
        return mol
    except Exception:
        return None


def load_pubchem_conformer(name):
    """Main helper: returns RDKit Mol of PubChem bioactive conformer if available"""
    cid = get_pubchem_cid(name)
    if cid is None:
        return None

    mol = load_pubchem_conformer_by_cid(cid)
    return mol
