# src/api/sync_chembl.py

from fastapi import APIRouter
from services.chembl_sync_service import sync_next_page, sync_targets_batch
from services.reactome_sync_service import sync_reactome_batch

router = APIRouter(prefix="/sync", tags=["chembl"])

@router.post("/reactome/pathways")
def sync_reactome_pathways(limit: int = 50):
    """
    Sync Reactome pathways for proteins that do not yet have pathways.

    - Uses UniProt IDs stored locally
    - Skips proteins already synced
    - Safe to call multiple times
    """

    result = sync_reactome_batch(limit=limit)
    return result


@router.post("/chembl/all")
def sync_chembl_all(
    max_drug_steps: int = 200,
    target_batch_size: int = 50
):
    """
    Orchestrates full ChEMBL sync:
    1. Approved drugs
    2. Targets & proteins
    """

    # ---- step 1: drugs ----
    final_drug_result = None
    for _ in range(max_drug_steps):
        result = sync_next_page()
        final_drug_result = result

        if result.get("status") == "completed":
            break
        if result.get("status") != "ok":
            break

    # ---- step 2: targets ----
    target_result = sync_targets_batch(limit=target_batch_size)

    return {
        "drugs": final_drug_result,
        "targets": target_result
    }


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