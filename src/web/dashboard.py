from fastapi import APIRouter, Request, Depends
from fastapi.templating import Jinja2Templates
from src.services.auth_dependencies import require_login

router = APIRouter()
templates = Jinja2Templates(directory="src/templates")


@router.get("/dashboard")
def dashboard(
    request: Request,
    user=Depends(require_login)
):
    return templates.TemplateResponse(
        "dashboard.html",
        {
            "request": request,
            "user": user,
        }
    )
