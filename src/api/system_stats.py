from fastapi import APIRouter, Depends
from src.services.auth_dependencies import require_login
from src.services.system_stats_service import (
    get_import_status_stats,
    get_not_in_chembl_reasons,
    get_local_chembl_count,
)

router = APIRouter(prefix="/system", tags=["system-stats"])


@router.get("/stats")
def system_stats(user=Depends(require_login)):
    """
    Returns system-wide statistics for dashboard / monitoring
    """

    return {
        "local_import_status": get_import_status_stats(),
        "not_in_chembl_reasons": get_not_in_chembl_reasons(),
        "chembl_local_count": get_local_chembl_count(),
    }
