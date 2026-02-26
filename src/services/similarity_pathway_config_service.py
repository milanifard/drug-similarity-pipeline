from sqlalchemy import text
from services.chemdb_service import get_engine


def get_selected_pathways():
    """
    Return active pathways configured for similarity table.
    """
    engine = get_engine()

    query = text("""
        SELECT pathway_id, display_order
        FROM similarity_selected_pathways
        WHERE is_active = TRUE
        ORDER BY display_order ASC
    """)

    with engine.connect() as conn:
        result = conn.execute(query)
        rows = result.fetchall()

    return [
        {
            "pathway_id": row[0],
            "display_order": row[1],
        }
        for row in rows
    ]


def add_selected_pathway(pathway_id: str, display_order: int = 0):
    engine = get_engine()

    query = text("""
        INSERT INTO similarity_selected_pathways (pathway_id, display_order, is_active)
        VALUES (:pathway_id, :display_order, TRUE)
    """)

    with engine.begin() as conn:
        conn.execute(query, {
            "pathway_id": pathway_id,
            "display_order": display_order
        })


def deactivate_selected_pathway(pathway_id: str):
    engine = get_engine()

    query = text("""
        UPDATE similarity_selected_pathways
        SET is_active = FALSE
        WHERE pathway_id = :pathway_id
    """)

    with engine.begin() as conn:
        conn.execute(query, {
            "pathway_id": pathway_id
        })