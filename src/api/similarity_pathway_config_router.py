from fastapi import APIRouter, Depends, Form, Request
from fastapi.templating import Jinja2Templates
from src.services.auth_dependencies import require_login

from services.similarity_pathway_config_service import (
    get_selected_pathways,
    add_selected_pathway,
    deactivate_selected_pathway,
)

router = APIRouter(
    prefix="/similarity-config",
    tags=["similarity-config"]
)

templates = Jinja2Templates(directory="src/templates")


# ✅ UI PAGE
@router.get("")
def similarity_config_page(
    request: Request,
    user=Depends(require_login)
):
    return templates.TemplateResponse(
        "similarity_config.html",
        {
            "request": request,
            "user": user
        }
    )


# ✅ API: list
@router.get("/pathways")
def list_selected_pathways(user=Depends(require_login)):
    pathways = get_selected_pathways()
    return {
        "count": len(pathways),
        "pathways": pathways
    }


# ✅ API: add
@router.post("/pathways")
def create_selected_pathway(
    pathway_id: str = Form(...),
    display_order: int = Form(0),
    user=Depends(require_login)
):
    add_selected_pathway(pathway_id, display_order)
    return {"status": "ok"}


# ✅ API: deactivate
@router.delete("/pathways/{pathway_id}")
def remove_selected_pathway(
    pathway_id: str,
    user=Depends(require_login)
):
    deactivate_selected_pathway(pathway_id)
    return {"status": "deactivated"}