from fastapi import APIRouter, Query, HTTPException, Depends
from src.services.similarity_service import compute_similarity_with_local_drugs
from src.services.auth_dependencies import require_login

router = APIRouter(prefix="/similarity", tags=["similarity"])


@router.get("/")
def similarity_with_local_drugs(
    drug: str = Query(...),
    alpha: float = Query(0.7, ge=0.0, le=1.0),
    radius: int = Query(3),
    nbits: int = Query(4096),
    threshold: float = Query(0.7, ge=0.0, le=1.0),
    user=Depends(require_login)
):
    try:
        df = compute_similarity_with_local_drugs(
            drug=drug,
            alpha=alpha,
            radius=radius,
            nbits=nbits,
        )
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except RuntimeError as e:
        raise HTTPException(status_code=500, detail=str(e))
    filtered = df[df["weighted_similarity"] >= threshold].reset_index(drop=True)
    return {
        "query": drug,
        "alpha": alpha,
        "threshold": threshold,
        "count": len(filtered),
        "results": filtered.to_dict(orient="records"),
    }
