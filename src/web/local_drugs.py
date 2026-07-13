from fastapi import APIRouter, Depends, Query, Request
from fastapi.templating import Jinja2Templates

from src.services.auth_dependencies import require_login
from src.services.chemdb_service import get_local_drugs_paginated


templates = Jinja2Templates(directory="src/templates")

router = APIRouter()


@router.get("/local-drugs")
def local_drugs_page(
    request: Request,
    user=Depends(require_login),
):
    """
    Display the local drugs list page.
    """

    return templates.TemplateResponse(
        "local_drugs.html",
        {
            "request": request,
            "user": user,
        },
    )


@router.get("/api/local-drugs")
def local_drugs_api(
    page: int = Query(default=1, ge=1),
    page_size: int = Query(default=20, ge=1, le=100),
    search: str = Query(default="", max_length=200),
    user=Depends(require_login),
):
    """
    Return paginated local drugs.

    Search is applied only to local_name.
    """

    return get_local_drugs_paginated(
        page=page,
        page_size=page_size,
        search=search,
    )