from fastapi import APIRouter, Request, Depends
from fastapi.templating import Jinja2Templates
from src.services.auth_dependencies import require_login

templates = Jinja2Templates(directory="src/templates")

router = APIRouter()


@router.get("/sync-chembl")
def sync_chembl_page(
    request: Request,
    user=Depends(require_login)
):
    return templates.TemplateResponse(
        "sync_chembl.html",
        {
            "request": request,
            "user": user,
        }
    )
