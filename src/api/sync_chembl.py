# src/api/sync_chembl.py

from fastapi import APIRouter
from services.chembl_sync_service import sync_next_page

router = APIRouter(prefix="/sync", tags=["chembl"])


@router.post("/chembl")
def sync_chembl(max_steps: int = 200):
    """
    Automatically loops sync process up to `max_steps` times,
    stopping when synchronization completes or an error occurs.
    """

    final_result = None

    for _ in range(max_steps):
        result = sync_next_page()
        final_result = result

        # stop if completed
        if result.get("status") == "completed":
            break

        # stop if something unexpected happened
        if result.get("status") != "ok":
            break

    return final_result