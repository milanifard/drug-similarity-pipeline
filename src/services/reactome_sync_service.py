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
            LEFT JOIN protein_reactome_status prs
            ON p.uniprot_id = prs.uniprot_id
            WHERE prs.uniprot_id IS NULL
            LIMIT :limit
        """), {"limit": limit}).fetchall()

    return [r[0] for r in rows]

def fetch_reactome_pathways(uniprot_id: str) -> List[Dict[str, Any]]:
    """
    Fetch Reactome pathways for a given UniProt ID.
    Raises exception on network/server errors.
    Returns empty list only if protein truly has no pathways.
    """
    url = f"{REACTOME_BASE_URL}/data/mapping/UniProt/{uniprot_id}/pathways"

    resp = requests.get(url, headers=HEADERS, timeout=TIMEOUT)

    # ---- VALID NO DATA CASE ----
    if resp.status_code == 404:
        return []

    # ---- RATE LIMIT ----
    if resp.status_code == 429:
        raise requests.exceptions.RequestException(
            f"Rate limited (429) for {uniprot_id}"
        )

    # ---- SERVER ERROR ----
    if 500 <= resp.status_code < 600:
        raise requests.exceptions.RequestException(
            f"Server error {resp.status_code} for {uniprot_id}"
        )

    # ---- OTHER BAD STATUS ----
    if resp.status_code != 200:
        raise requests.exceptions.RequestException(
            f"Unexpected status {resp.status_code} for {uniprot_id}"
        )

    return resp.json()


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

def save_reactome_status(uniprot_id: str, has_pathway: bool):
    engine = get_engine()
    with engine.begin() as conn:
        conn.execute(text("""
            INSERT IGNORE INTO protein_reactome_status
                (uniprot_id, has_pathway)
            VALUES
                (:uid, :has_pathway)
        """), {
            "uid": uniprot_id,
            "has_pathway": has_pathway
        })

def sync_protein_pathways(uniprot_id: str):
    try:
        pathways = fetch_reactome_pathways(uniprot_id)
    except Exception as e:
            print(f"[FETCH ERROR] {uniprot_id}: {e}")
            return
    if pathways is None:
        return

    if len(pathways) == 0:
        save_reactome_status(uniprot_id, has_pathway=False)
        return
    
    save_pathways(pathways)

    pathway_ids = [
        p.get("stId")
        for p in pathways
        if p.get("stId")
    ]

    save_protein_pathways(uniprot_id, pathway_ids)
    save_reactome_status(uniprot_id, has_pathway=True)

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
