from sqlalchemy import text
from services.chemdb_service import get_engine


def get_selected_targets():
    engine = get_engine()
    query = text("""
        SELECT target_chembl_id, target_name, display_order
        FROM similarity_selected_targets
        WHERE is_active = TRUE
        ORDER BY display_order ASC, target_name ASC
    """)

    with engine.connect() as conn:
        rows = conn.execute(query).fetchall()

    return [
        {
            "target_chembl_id": row[0],
            "target_name": row[1],
            "display_order": row[2],
        }
        for row in rows
    ]


def add_selected_target(
    target_chembl_id: str,
    target_name: str | None = None,
    display_order: int = 0,
):
    engine = get_engine()

    query = text("""
        INSERT INTO similarity_selected_targets
            (target_chembl_id, target_name, display_order, is_active)
        VALUES
            (:target_chembl_id, :target_name, :display_order, TRUE)
        ON DUPLICATE KEY UPDATE
            target_name = VALUES(target_name),
            display_order = VALUES(display_order),
            is_active = TRUE
    """)

    with engine.begin() as conn:
        conn.execute(query, {
            "target_chembl_id": target_chembl_id,
            "target_name": target_name,
            "display_order": display_order,
        })


def deactivate_selected_target(target_chembl_id: str):
    engine = get_engine()

    query = text("""
        UPDATE similarity_selected_targets
        SET is_active = FALSE
        WHERE target_chembl_id = :target_chembl_id
    """)

    with engine.begin() as conn:
        conn.execute(query, {"target_chembl_id": target_chembl_id})