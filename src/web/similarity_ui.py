from fastapi import APIRouter, Request, Form, Depends
from fastapi.templating import Jinja2Templates
from src.services.auth_dependencies import require_login
from src.services.similarity_service import compute_similarity_with_local_drugs

router = APIRouter()
templates = Jinja2Templates(directory="src/templates")


@router.get("/similarity-ui")
def similarity_page(
    request: Request,
    user=Depends(require_login)
):
    return templates.TemplateResponse(
        "similarity.html",
        {
            "request": request,
            "user": user,   # ✅ خیلی مهم
        }
    )


@router.post("/similarity-ui/run")
def run_similarity(
    request: Request,
    user=Depends(require_login),   # حتی اگر استفاده نشود
    drug: str = Form(...),
    alpha: float = Form(0.7),
    threshold: float = Form(0.5),
    radius: int = Form(3),
    nbits: int = Form(4096),
):
    df = compute_similarity_with_local_drugs(
        drug=drug,
        alpha=alpha,
        radius=radius,
        nbits=nbits,
    )

    df = df[df["weighted_similarity"] >= threshold]
    df = df.round(3)

    return {
        "count": len(df),
        "results": df.to_dict(orient="records"),
    }
