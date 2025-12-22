from fastapi import APIRouter, HTTPException
from services.drug_pathway_service import get_drug_pathways

router = APIRouter(prefix="/drugs", tags=["drug-pathways"])


@router.get("/{chembl_id}/pathways")
def drug_pathways(chembl_id: str):
    """
    Get Reactome pathways associated with a ChEMBL drug ID.
    """

    pathways = get_drug_pathways(chembl_id)

    if not pathways:
        raise HTTPException(
            status_code=404,
            detail=f"No pathways found for drug {chembl_id}"
        )

    return {
        "chembl_id": chembl_id,
        "pathway_count": len(pathways),
        "pathways": pathways
    }
