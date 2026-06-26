from fastapi import APIRouter, Depends, Form

from src.services.auth_dependencies import require_login
from services.similarity_target_config_service import (
    get_selected_targets,
    add_selected_target,
    deactivate_selected_target,
)

router = APIRouter(
    prefix="/similarity-config",
    tags=["similarity-target-config"],
)


@router.get("/targets")
def list_selected_targets(user=Depends(require_login)):
    targets = get_selected_targets()
    return {
        "count": len(targets),
        "targets": targets,
    }


@router.post("/targets")
def create_selected_target(
    target_chembl_id: str = Form(...),
    target_name: str = Form(None),
    display_order: int = Form(0),
    user=Depends(require_login),
):
    add_selected_target(
        target_chembl_id=target_chembl_id,
        target_name=target_name,
        display_order=display_order,
    )
    return {"status": "ok"}


@router.delete("/targets/{target_chembl_id}")
def remove_selected_target(
    target_chembl_id: str,
    user=Depends(require_login),
):
    deactivate_selected_target(target_chembl_id)
    return {"status": "deactivated"}