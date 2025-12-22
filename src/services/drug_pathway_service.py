from typing import List, Dict, Any
from sqlalchemy import text
from services.chemdb_service import get_engine


def get_drug_pathways(molecule_chembl_id: str) -> List[Dict[str, Any]]:
    """
    Return Reactome pathways associated with a given ChEMBL drug ID
    using locally cached data only.
    """

    engine = get_engine()

    sql = text("""
        SELECT DISTINCT
            rp.pathway_id,
            rp.pathway_name,
            rp.top_level_class,
            rp.species
        FROM drug_targets dt
        JOIN target_proteins tp
            ON dt.target_chembl_id = tp.target_chembl_id
        JOIN protein_pathways pp
            ON tp.uniprot_id = pp.uniprot_id
        JOIN reactome_pathways rp
            ON pp.pathway_id = rp.pathway_id
        WHERE dt.molecule_chembl_id = :chembl_id
        ORDER BY rp.pathway_name
    """)

    with engine.begin() as conn:
        rows = conn.execute(sql, {"chembl_id": molecule_chembl_id}).mappings().all()

    return [dict(r) for r in rows]
