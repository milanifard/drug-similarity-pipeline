from fastapi import APIRouter, Request, UploadFile, File, Form, Depends
from fastapi.templating import Jinja2Templates
from src.services.auth_dependencies import require_login
from src.api.import_local_drugs import import_local_drugs

router = APIRouter()
templates = Jinja2Templates(directory="src/templates")


@router.get("/import")
def import_page(
    request: Request,
    user=Depends(require_login)
):
    return templates.TemplateResponse(
        "import.html",
        {"request": request, "user": user}
    )
