# chemdb_service.py

import os
from typing import Dict, Any, Iterable
from sqlalchemy import create_engine, text
from sqlalchemy.engine import Engine


def get_engine() -> Engine:
    db_url = os.getenv(
        "CHEMDB_URL",
        "mysql+pymysql://chemuser:chempass@chemdb:3306/chemdb"
    )
    return create_engine(db_url, pool_pre_ping=True, future=True)


# ----------------------------------------------------
# local_drugs table CRUD
# ----------------------------------------------------
def bulk_upsert_local_drugs(rows: Iterable[Dict[str, Any]]):
    """
    Insert or update rows inside local_drugs table.
    """
    sql = text("""
        INSERT INTO local_drugs
            (normalized_name, local_name, chembl_name, molecule_chembl_id,
             smiles, conformer_source, molblock)
        VALUES
            (:normalized_name, :local_name, :chembl_name, :molecule_chembl_id,
             :smiles, :conformer_source, :molblock)
        ON DUPLICATE KEY UPDATE
            local_name = VALUES(local_name),
            chembl_name = VALUES(chembl_name),
            molecule_chembl_id = VALUES(molecule_chembl_id),
            smiles = VALUES(smiles),
            conformer_source = VALUES(conformer_source),
            molblock = VALUES(molblock);
    """)

    engine = get_engine()
    with engine.begin() as conn:
        conn.execute(sql, list(rows))


def get_local_drug_by_normalized_name(normalized_name: str):
    """
    Fetch a single stored molecule by normalized name.
    """
    sql = text("""
        SELECT id, normalized_name, local_name, chembl_name,
               molecule_chembl_id, smiles, conformer_source, molblock
        FROM local_drugs
        WHERE normalized_name = :name
        LIMIT 1
    """)

    engine = get_engine()
    with engine.connect() as conn:
        row = conn.execute(sql, {"name": normalized_name}).mappings().first()
        return dict(row) if row else None

def get_local_drugs_paginated(
    page: int = 1,
    page_size: int = 20,
    search: str = "",
) -> Dict[str, Any]:
    """
    Return paginated local drugs.

    Search is applied only to local_name.
    """

    page = max(page, 1)
    page_size = min(max(page_size, 1), 100)
    offset = (page - 1) * page_size

    search = search.strip()

    where_clause = ""
    params: Dict[str, Any] = {
        "limit": page_size,
        "offset": offset,
    }

    if search:
        where_clause = "WHERE local_name LIKE :search"
        params["search"] = f"%{search}%"

    count_sql = text(f"""
        SELECT COUNT(*) AS total
        FROM local_drugs
        {where_clause}
    """)

    list_sql = text(f"""
        SELECT
            local_name,
            chembl_name,
            molecule_chembl_id
        FROM local_drugs
        {where_clause}
        ORDER BY local_name ASC
        LIMIT :limit OFFSET :offset
    """)

    engine = get_engine()

    with engine.connect() as conn:
        total = conn.execute(count_sql, params).scalar_one()

        rows = conn.execute(
            list_sql,
            params,
        ).mappings().all()

    total_pages = (
        (total + page_size - 1) // page_size
        if total > 0
        else 0
    )

    return {
        "items": [dict(row) for row in rows],
        "page": page,
        "page_size": page_size,
        "total": total,
        "total_pages": total_pages,
        "search": search,
    }