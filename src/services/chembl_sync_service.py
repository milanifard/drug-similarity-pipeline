# src/service/chembl_sync_service.py

from typing import List, Dict, Any
from sqlalchemy import text
from services.chemdb_service import get_engine
from data_sources.chembl_client import fetch_approved_page


PAGE_SIZE = 20   # ChEMBL returns ~20 records per page


def get_sync_state() -> Dict[str, Any]:
    """Load current sync state from DB."""
    engine = get_engine()
    with engine.begin() as conn:
        row = conn.execute(text("""
            SELECT last_offset, is_completed
            FROM chembl_sync_state
            WHERE id = 1
        """)).mappings().first()

        if row:
            return dict(row)

        # create initial state
        conn.execute(text("""
            INSERT INTO chembl_sync_state (id, last_offset, is_completed)
            VALUES (1, 0, 0)
        """))

        return {"last_offset": 0, "is_completed": False}


def update_sync_state(offset: int, completed: bool):
    """Update sync progress."""
    engine = get_engine()
    with engine.begin() as conn:
        result = conn.execute(text("""
            UPDATE chembl_sync_state
            SET last_offset = :off,
                is_completed = :done
            WHERE id = 1
        """), {"off": offset, "done": completed})

        # if row does not exist → insert
        if result.rowcount == 0:
            conn.execute(text("""
                INSERT INTO chembl_sync_state (id, last_offset, is_completed)
                VALUES (1, :off, :done)
            """), {"off": offset, "done": completed})


def save_chembl_page(rows: List[Dict[str, Any]]):
    """Insert or ignore (if duplicate)."""
    if not rows:
        return

    engine = get_engine()

    sql = text("""
        INSERT IGNORE INTO chembl_approved_drugs
            (molecule_chembl_id, name, smiles)
        VALUES
            (:chembl_id, :name, :smiles)
    """)

    with engine.begin() as conn:
        conn.execute(sql, rows)


def sync_next_page() -> Dict[str, Any]:
    """
    Fetch the next page from ChEMBL and store locally.
    Returns status info for the API.
    """

    state = get_sync_state()

    if state["is_completed"]:
        return {
            "status": "completed",
            "message": "All ChEMBL approved drugs already synchronized."
        }

    offset = state["last_offset"]

    print(f"[SYNC] Fetching ChEMBL page offset={offset}")

    # Fetch next 20 records
    page = fetch_approved_page(offset=offset, limit=PAGE_SIZE)

    if len(page) == 0:
        update_sync_state(offset, True)
        return {
            "status": "completed",
            "message": f"Sync finished at offset={offset}"
        }

    # Save in DB
    save_chembl_page(page)

    # Update state
    update_sync_state(offset + PAGE_SIZE, False)

    return {
        "status": "ok",
        "fetched": len(page),
        "new_offset": offset + PAGE_SIZE
    }
