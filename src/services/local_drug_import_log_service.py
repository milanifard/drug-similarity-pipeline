from sqlalchemy import text
from src.services.chemdb_service import get_engine


def log_import_event(
    normalized_name: str,
    status: str,
    raw_name: str | None = None,
    chembl_id: str | None = None,
    chembl_name: str | None = None,
    smiles: str | None = None,
    num_atoms: int | None = None,
    message: str | None = None,
):
    """
    Insert a log record only if (normalized_name, status) does not already exist.
    """
    engine = get_engine()

    query = text("""
        INSERT INTO local_drug_import_log
        (normalized_name, raw_name, chembl_id, chembl_name, smiles,
         status, message, num_atoms)
        VALUES
        (:normalized_name, :raw_name, :chembl_id, :chembl_name, :smiles,
         :status, :message, :num_atoms)
        ON DUPLICATE KEY UPDATE
            message = VALUES(message)
    """)

    with engine.begin() as conn:
        conn.execute(query, {
            "normalized_name": normalized_name,
            "raw_name": raw_name,
            "chembl_id": chembl_id,
            "chembl_name": chembl_name,
            "smiles": smiles,
            "status": status,
            "message": message,
            "num_atoms": num_atoms,
        })
