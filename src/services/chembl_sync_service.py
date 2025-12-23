# src/service/chembl_sync_service.py

from typing import List, Dict, Any, Set
from sqlalchemy import text
from services.chemdb_service import get_engine
from data_sources.chembl_client import fetch_approved_page, fetch_target_components, fetch_targets_for_drug


PAGE_SIZE = 20   # ChEMBL returns ~20 records per page

def sync_targets_batch(limit: int = 50) -> Dict[str, Any]:
    """
    Sync targets for a batch of local drugs.
    """
    drug_ids = get_unsynced_drug_ids(limit=limit)

    if not drug_ids:
        return {
            "status": "completed",
            "message": "No drugs found to sync targets."
        }

    for chembl_id in drug_ids:
        print(f"[TARGET-SYNC] {chembl_id}")
        try:
            sync_drug_targets(chembl_id)
        except Exception as ex:
            print(f"[TARGET-SYNC ERROR] {chembl_id}: {ex}")

    return {
        "status": "ok",
        "processed_drugs": len(drug_ids)
    }

def get_unsynced_drug_ids(limit: int = 50) -> List[str]:
    """
    Return drugs whose targets are not yet synced.
    """
    engine = get_engine()
    with engine.begin() as conn:
        rows = conn.execute(
            text("""
                SELECT molecule_chembl_id
                FROM local_drugs
                WHERE targets_synced = 0
                LIMIT :limit
            """),
            {"limit": limit}
        ).fetchall()

    return [r[0] for r in rows]


def drug_targets_exist(molecule_chembl_id: str) -> bool:
    engine = get_engine()
    with engine.begin() as conn:
        row = conn.execute(text("""
            SELECT 1
            FROM drug_targets
            WHERE molecule_chembl_id = :cid
            LIMIT 1
        """), {"cid": molecule_chembl_id}).fetchone()
    return row is not None

def target_proteins_exist(target_chembl_id: str) -> bool:
    engine = get_engine()
    with engine.begin() as conn:
        row = conn.execute(text("""
            SELECT 1
            FROM target_proteins
            WHERE target_chembl_id = :tid
            LIMIT 1
        """), {"tid": target_chembl_id}).fetchone()
    return row is not None


def get_local_drug_ids(limit: int = 100) -> List[str]:
    """
    Get a batch of molecule_chembl_id that are already synced locally.
    """
    engine = get_engine()
    with engine.begin() as conn:
        rows = conn.execute(text("""
            SELECT molecule_chembl_id
            FROM chembl_approved_drugs
            ORDER BY molecule_chembl_id
            LIMIT :limit
        """), {"limit": limit}).fetchall()

    return [r[0] for r in rows]

def save_targets(targets: List[Dict[str, Any]]):
    """
    Insert targets into chembl_targets table.
    """
    if not targets:
        return

    engine = get_engine()
    sql = text("""
        INSERT IGNORE INTO chembl_targets
            (target_chembl_id, target_name, organism, target_type)
        VALUES
            (:target_chembl_id, :target_name, :organism, :target_type)
    """)

    with engine.begin() as conn:
        conn.execute(sql, targets)

def save_drug_targets(molecule_chembl_id: str, target_ids: Set[str]):
    """
    Save mapping between drug and targets.
    """
    if not target_ids:
        return

    rows = [
        {
            "molecule_chembl_id": molecule_chembl_id,
            "target_chembl_id": tid
        }
        for tid in target_ids
    ]

    engine = get_engine()
    sql = text("""
        INSERT IGNORE INTO drug_targets
            (molecule_chembl_id, target_chembl_id)
        VALUES
            (:molecule_chembl_id, :target_chembl_id)
    """)

    with engine.begin() as conn:
        conn.execute(sql, rows)

def save_proteins(proteins: List[Dict[str, Any]]):
    """
    Insert proteins into proteins table.
    """
    if not proteins:
        return

    engine = get_engine()
    sql = text("""
        INSERT IGNORE INTO proteins
            (uniprot_id, protein_name, gene_name, organism)
        VALUES
            (:uniprot_id, :protein_name, :gene_name, :organism)
    """)

    with engine.begin() as conn:
        conn.execute(sql, proteins)

def save_target_proteins(target_chembl_id: str, uniprot_ids: Set[str]):
    """
    Save mapping between target and proteins.
    """
    if not uniprot_ids:
        return

    rows = [
        {
            "target_chembl_id": target_chembl_id,
            "uniprot_id": uid
        }
        for uid in uniprot_ids
    ]

    engine = get_engine()
    sql = text("""
        INSERT IGNORE INTO target_proteins
            (target_chembl_id, uniprot_id)
        VALUES
            (:target_chembl_id, :uniprot_id)
    """)

    with engine.begin() as conn:
        conn.execute(sql, rows)

def mark_drug_targets_synced(chembl_id: str):
    engine = get_engine()
    with engine.begin() as conn:
        conn.execute(
            text("""
                UPDATE local_drugs
                SET targets_synced = 1
                WHERE molecule_chembl_id = :id
            """),
            {"id": chembl_id}
        )

def sync_drug_targets(molecule_chembl_id: str):
    """
    Sync targets and proteins for a single drug.
    """
    if drug_targets_exist(molecule_chembl_id):
            print(f"[SKIP] targets already synced for {molecule_chembl_id}")
            return
    
    targets = fetch_targets_for_drug(molecule_chembl_id)

    if not targets:
        return

    # ---- save targets ----
    save_targets(targets)

    target_ids = {t["target_chembl_id"] for t in targets}
    save_drug_targets(molecule_chembl_id, target_ids)
    mark_drug_targets_synced(molecule_chembl_id)
    # ---- for each target → proteins ----
    for t in targets:
        target_id = t["target_chembl_id"]

        components = fetch_target_components(target_id)
        if not components:
            continue

        save_proteins(components)

        uniprot_ids = {
            c["uniprot_id"]
            for c in components
            if c.get("uniprot_id")
        }

        save_target_proteins(target_id, uniprot_ids)



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
