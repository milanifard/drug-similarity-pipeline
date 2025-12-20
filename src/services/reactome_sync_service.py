import requests
from typing import List, Dict, Any
from sqlalchemy import text
from services.chemdb_service import get_engine

REACTOME_BASE_URL = "https://reactome.org/ContentService"
HEADERS = {
    "Accept": "application/json"
}
TIMEOUT = 20

def get_unsynced_uniprot_ids(limit: int = 100) -> List[str]:
    """
    Return UniProt IDs that exist locally but have no pathways yet.
    """
    engine = get_engine()
    with engine.begin() as conn:
        rows = conn.execute(text("""
            SELECT p.uniprot_id
            FROM proteins p
            LEFT JOIN protein_pathways pp
              ON p.uniprot_id = pp.uniprot_id
            WHERE pp.uniprot_id IS NULL
            LIMIT :limit
        """), {"limit": limit}).fetchall()

    return [r[0] for r in rows]

def fetch_reactome_pathways(uniprot_id: str) -> List[Dict[str, Any]]:
    """
    Fetch Reactome pathways for a given UniProt ID.
    """
    url = f"{REACTOME_BASE_URL}/data/mapping/UniProt/{uniprot_id}/pathways"

    try:
        resp = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
        if resp.status_code != 200:
            print(f"[REACTOME ERROR] {uniprot_id} → HTTP {resp.status_code}")
            return []

        return resp.json()

    except Exception as ex:
        print(f"[REACTOME ERROR] {uniprot_id}: {ex}")
        return []


def save_pathways(pathways: List[Dict[str, Any]]):
    if not pathways:
        return

    rows = []
    for p in pathways:
        rows.append({
            "pathway_id": p.get("stId"),
            "pathway_name": p.get("displayName"),
            "top_level_class": p.get("topLevelPathway", {}).get("displayName")
                                if p.get("topLevelPathway") else None,
            "species": p.get("speciesName"),
        })

    engine = get_engine()
    sql = text("""
        INSERT IGNORE INTO reactome_pathways
            (pathway_id, pathway_name, top_level_class, species)
        VALUES
            (:pathway_id, :pathway_name, :top_level_class, :species)
    """)

    with engine.begin() as conn:
        conn.execute(sql, rows)

def save_protein_pathways(uniprot_id: str, pathway_ids: List[str]):
    if not pathway_ids:
        return

    rows = [
        {
            "uniprot_id": uniprot_id,
            "pathway_id": pid
        }
        for pid in pathway_ids
    ]

    engine = get_engine()
    sql = text("""
        INSERT IGNORE INTO protein_pathways
            (uniprot_id, pathway_id)
        VALUES
            (:uniprot_id, :pathway_id)
    """)

    with engine.begin() as conn:
        conn.execute(sql, rows)

def sync_protein_pathways(uniprot_id: str):
    pathways = fetch_reactome_pathways(uniprot_id)

    if not pathways:
        return

    save_pathways(pathways)

    pathway_ids = [
        p.get("stId")
        for p in pathways
        if p.get("stId")
    ]

    save_protein_pathways(uniprot_id, pathway_ids)

def sync_reactome_batch(limit: int = 50) -> Dict[str, Any]:
    """
    Sync Reactome pathways for proteins not yet processed.
    """
    uniprot_ids = get_unsynced_uniprot_ids(limit=limit)

    if not uniprot_ids:
        return {
            "status": "completed",
            "message": "All proteins already have Reactome pathways."
        }

    for uid in uniprot_ids:
        print(f"[REACTOME-SYNC] {uid}")
        sync_protein_pathways(uid)

    return {
        "status": "ok",
        "processed_proteins": len(uniprot_ids)
    }
