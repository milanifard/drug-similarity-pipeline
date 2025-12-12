# src/services/drug_name_normalizer.py

import re
from typing import Dict, List, Set

# Allowed short abbreviations
SHORT_WHITELIST: Set[str] = {"bcg", "hcg", "fsh", "lh", "tsh", "gcsf"}

# Banned non-drug / formulation terms (conservative)
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


def normalize_name(raw: str) -> List[str]:
    """
    Low-level normalization for a single raw drug name:
      - Remove parentheses
      - Split combination drugs by '+'
      - Convert to lowercase
      - Remove common salts (hydrochloride, sulfate, ...)
      - Remove dosage forms (tablet, injection, ...)
      - Keep only letters/numbers/space/hyphen
    Returns a list of normalized tokens (for combinations).
    """
    if not isinstance(raw, str):
        return []

    # Remove parenthesized parts
    name = re.sub(r"\(.*?\)", "", raw)

    # Split combination drugs by '+'
    parts = re.split(r"\s*\+\s*", name)

    cleaned: List[str] = []

    for p in parts:
        p = p.strip().lower()
        if not p:
            continue

        # Remove salt descriptors
        p = re.sub(
            r"\b(as|as\.)?\s*(hydrochloride|hcl|sulfate|sulphate|na|sodium|"
            r"potassium|phosphate|carbonate|mesilate|mesylate|acetate|tartrate|"
            r"nitrate|maleate|fumarate|succinate|bromide|iodide)\b",
            "",
            p,
        )

        # Remove dosage forms
        p = re.sub(
            r"\b(tablet|tablets|capsule|capsules|injection|injections|solution|"
            r"suspension|cream|ointment|gel|drops|syrup|patch|spray)\b",
            "",
            p,
        )

        # Keep only [a-z0-9 -]
        p = re.sub(r"[^a-z0-9\s\-]", "", p)
        p = re.sub(r"\s+", " ", p).strip()

        if p:
            cleaned.append(p)

    return cleaned


def is_radiopharmaceutical(name: str) -> bool:
    """
    Detect radiopharmaceutical-like patterns:
      - Starting with digits + letters + '-'
      - Patterns like [131i], [32p], ...
    """
    if re.match(r"^\d+\s*[a-z]*-", name):
        return True
    if re.search(r"\[\d+\s*[a-z]+\]", name):
        return True
    return False


def is_banned_term(name: str) -> bool:
    """
    Return True if a normalized name clearly belongs to
    non-small-molecule categories (vaccines, blood products, etc.).
    """
    for key in BANNED_KEYWORDS:
        if key in name:
            return True
    return False


def normalize_raw_names(raw_names: List[str]) -> List[str]:
    """
    High-level normalization for a list of raw drug names (conservative).

    Filters:
      - length <= 3 (unless whitelisted)
      - radiopharmaceutical-like patterns
      - banned non-drug categories
    Returns a sorted unique list of normalized names.
    """
    normalized: List[str] = []

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


def build_normalized_index(names: List[str]) -> Dict[str, List[str]]:
    """
    Build mapping:
        normalized_name -> list of original names

    Uses the same normalization rules as normalize_raw_names()
    (except banning, which usually is applied on the caller side).
    """
    index: Dict[str, List[str]] = {}

    for raw in names:
        for norm in normalize_name(raw):
            if is_radiopharmaceutical(norm):
                continue

            if len(norm) <= 3 and norm not in SHORT_WHITELIST:
                continue

            index.setdefault(norm, []).append(raw)

    return index
