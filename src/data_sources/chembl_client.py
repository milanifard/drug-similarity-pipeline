from chembl_webresource_client.new_client import new_client as _new_client
import pandas as pd
from sqlalchemy import text
from services.chemdb_service import get_engine
from chembl_webresource_client.http_errors import HttpApplicationError

_client = None

def get_client():
    global _client
    if _client is None:
        _client = _new_client  # initialize only on first use
    return _client

def get_targets(disease_name):
    """Search biological targets associated with a given disease name."""
    client = get_client()
    target = client.target
    targets = target.search(disease_name)
    return pd.DataFrame(targets)


def get_activities(target_chembl_id):
    """Fetch activity records associated with a specific biological target."""
    client = get_client()
    activity = client.activity
    res = activity.filter(target_chembl_id=target_chembl_id).only(
        ['molecule_chembl_id', 'pchembl_value']
    )
    return pd.DataFrame(res)


def get_smiles(molecule_ids, limit=10):
    """Fetch SMILES strings for a list of ChEMBL molecule IDs."""
    client = get_client()
    molecule = client.molecule
    smiles_list = []
    for mol_id in molecule_ids[:limit]:
        mol = molecule.get(mol_id)
        if mol and mol['molecule_structures']:
            smiles_list.append(mol['molecule_structures']['canonical_smiles'])
    return smiles_list


def search_molecule(drug_name: str):
    """Search molecular information based on a drug name (e.g., 'Memantine')."""
    client = get_client()
    molecule = client.molecule
    res = molecule.filter(pref_name__icontains=drug_name).only(
        ['molecule_chembl_id', 'pref_name', 'max_phase', 'molecule_structures']
    )

    data = []
    for r in res:
        if r.get('molecule_structures'):
            data.append({
                'chembl_id': r['molecule_chembl_id'],
                'name': r.get('pref_name', ''),
                'max_phase': r.get('max_phase', ''),
                'smiles': r['molecule_structures']['canonical_smiles']
            })
    return pd.DataFrame(data)


def fetch_approved_page(offset: int, limit: int = 20):
    """
    Fetch a single page of phase-4 (approved) molecules.
    Returns a list of dicts: [{chembl_id, name, smiles}, ...]
    """
    client = get_client()
    molecule = client.molecule

    try:
        page = molecule.filter(max_phase=4).only(
            ['molecule_chembl_id', 'pref_name', 'molecule_structures']
        )[offset: offset + limit]
    except Exception as ex:
        print(f"[CHEMBL ERROR] Paging failed at offset {offset}: {ex}")
        return []

    rows = []
    for r in page:
        if r.get("molecule_structures"):
            rows.append({
                "chembl_id": r["molecule_chembl_id"],
                "name": r.get("pref_name", ""),
                "smiles": r["molecule_structures"]["canonical_smiles"],
            })

    return rows

def get_local_approved_drugs(limit: int | None = None) -> pd.DataFrame:
    """
    Return approved ChEMBL drugs stored in the local cache table.
    Does NOT fetch anything from ChEMBL.
    """

    engine = get_engine()

    query = "SELECT molecule_chembl_id, name, smiles, created_at FROM chembl_approved_drugs ORDER BY id"

    if limit is not None:
        query += f" LIMIT {limit}"

    return pd.read_sql(query, engine)
