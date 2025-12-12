from typing import List
from services.local_drug_service import SHORT_WHITELIST, normalize_name, is_radiopharmaceutical

def normalize_raw_names(raw_names: List[str]) -> List[str]:
    """
    Normalize a list of raw drug names and return a unique clean list.
    Filters:
      - length <= 3 (unless whitelisted)
      - radiopharmaceutical patterns
    """
    normalized = []

    for raw in raw_names:
        for norm in normalize_name(raw):
            if len(norm) <= 3 and norm not in SHORT_WHITELIST:
                continue
            if is_radiopharmaceutical(norm):
                continue
            normalized.append(norm)

    return sorted(set(normalized))
