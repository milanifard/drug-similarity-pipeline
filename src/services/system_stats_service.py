from sqlalchemy import text
from src.services.chemdb_service import get_engine


def get_import_status_stats():
    """
    Returns count of local drug import records grouped by status
    """
    engine = get_engine()

    sql = text("""
        SELECT status, COUNT(*) AS count
        FROM local_drug_import_log
        GROUP BY status
    """)

    with engine.connect() as conn:
        rows = conn.execute(sql).fetchall()

    return [
        {"status": row.status, "count": row.count}
        for row in rows
    ]


def get_not_in_chembl_reasons():
    """
    Returns reasons (message) for drugs not found in ChEMBL with counts
    """
    engine = get_engine()

    sql = text("""
        SELECT message, COUNT(*) AS count
        FROM local_drug_import_log
        WHERE status = 'NOT_IN_CHEMBL'
        GROUP BY message
        ORDER BY count DESC
    """)

    with engine.connect() as conn:
        rows = conn.execute(sql).fetchall()

    return [
        {"message": row.message, "count": row.count}
        for row in rows
    ]


def get_local_chembl_count():
    """
    Returns total number of locally cached approved ChEMBL drugs
    """
    engine = get_engine()

    sql = text("SELECT COUNT(*) AS count FROM chembl_approved_drugs")

    with engine.connect() as conn:
        row = conn.execute(sql).first()

    return row.count if row else 0
