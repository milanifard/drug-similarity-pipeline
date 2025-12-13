import pandas as pd
from sqlalchemy import text
from services.chemdb_service import get_engine
from chembl_webresource_client.http_errors import HttpApplicationError

_client = None


def get_client():
    """
    Lazy-load ChEMBL client.
    Importing chembl_webresource_client at module-load time causes
    service startup failures if the API is down.
    """

    global _client
    if _client is None:
        try:
            # 🔥 Lazy import — ONLY happens when first needed
            from chembl_webresource_client.new_client import new_client as _new_client
            _client = _new_client
        except Exception as ex:
            raise RuntimeError(f"Failed to initialize ChEMBL client: {ex}")

    return _client


# ============================================================
# SAFE API CALL WRAPPERS
# ============================================================

def _safe_call(func, *args, **kwargs):
    """
    Universal wrapper for ChEMBL calls to avoid crashing the service.
    Returns empty DataFrame or empty list on failure.
    """
    try:
        return func(*args, **kwargs)
    except Exception as ex:
        print(f"[CHEMBL ERROR] {ex}")
        return None


# ============================================================
# API METHODS
# ============================================================

def get_targets(disease_name):
    client = get_client()
    target = client.target

    res = _safe_call(target.search, disease_name)
    return pd.DataFrame(res or [])


def get_activities(target_chembl_id):
    client = get_client()
    activity = client.activity

    try:
        res = activity.filter(target_chembl_id=target_chembl_id).only(
            ['molecule_chembl_id', 'pchembl_value']
        )
        return pd.DataFrame(res)
    except Exception as ex:
        print(f"[CHEMBL ERROR] {ex}")
        return pd.DataFrame()


def get_smiles(molecule_ids, limit=10):
    client = get_client()
    molecule = client.molecule

    smiles_list = []
    for mol_id in molecule_ids[:limit]:
        try:
            mol = molecule.get(mol_id)
            if mol and mol.get("molecule_structures"):
                smiles_list.append(mol["molecule_structures"]["canonical_smiles"])
        except Exception as ex:
            print(f"[CHEMBL ERROR] get_smiles({mol_id}): {ex}")
            continue

    return smiles_list


def search_molecule(drug_name: str):
    client = get_client()
    molecule = client.molecule

    try:
        res = molecule.filter(pref_name__icontains=drug_name).only(
            ['molecule_chembl_id', 'pref_name', 'max_phase', 'molecule_structures']
        )
    except Exception as ex:
        print(f"[CHEMBL ERROR] search_molecule({drug_name}): {ex}")
        return pd.DataFrame()

    data = []
    for r in res:
        if r.get("molecule_structures"):
            data.append({
                'chembl_id': r['molecule_chembl_id'],
                'name': r.get('pref_name', ''),
                'max_phase': r.get('max_phase', ''),
                'smiles': r['molecule_structures']['canonical_smiles']
            })

    return pd.DataFrame(data)


def fetch_approved_page(offset: int, limit: int = 20):
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


# ============================================================
# LOCAL CACHE METHODS
# ============================================================

def get_local_approved_drugs(limit: int | None = None) -> pd.DataFrame:
    engine = get_engine()

    query = "SELECT molecule_chembl_id, name, smiles, created_at FROM chembl_approved_drugs ORDER BY id"

    if limit is not None:
        query += f" LIMIT {limit}"

    return pd.read_sql(query, engine)
