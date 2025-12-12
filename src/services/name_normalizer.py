# src/services/name_normalizer.py

from typing import List
from services.local_drug_service import SHORT_WHITELIST, normalize_name, is_radiopharmaceutical


# --- NEW: banned patterns (non-drugs, vaccines, blood products, etc.) ---
BANNED_KEYWORDS = [
    "vaccine", "virus", "immune", "immunoglobulin", "globulin",
    "plasma", "platelet", "blood", "erythrocyte",
    "dialysis", "hemodialysis", "peritoneal",
    "electrolyte", "solution", "infusion",
    "nutrition", "parenteral", "multivitamin", "vitamin",
    "bath", "ointment", "oil",
    "serum", "antitoxin", "antivenom",
    "test", "skin test",
    "saline", "dextrose", "glucose", "nacl",
    "charcoal", "talc", "gelatin",
    "oxygen", "nitrous", "gas",
]


def is_banned_term(name: str) -> bool:
    """Return True if a normalized name contains any banned keyword."""
    for key in BANNED_KEYWORDS:
        if key in name:
            return True
    return False


def normalize_raw_names(raw_names: List[str]) -> List[str]:
    """
    Normalize a list of raw drug names (conservative).
    Filters:
      - length <= 3 (unless whitelisted)
      - radiopharmaceuticals
      - banned non-drug categories
    """
    normalized = []

    for raw in raw_names:
        for norm in normalize_name(raw):
            if len(norm) <= 3 and norm not in SHORT_WHITELIST:
                continue

            if is_radiopharmaceutical(norm):
                continue

            if is_banned_term(norm):
                continue

            normalized.append(norm)

    return sorted(set(normalized))
