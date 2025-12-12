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
