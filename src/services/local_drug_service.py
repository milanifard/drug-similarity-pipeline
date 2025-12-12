# local_drug_service.py

import io
import re
from typing import Dict, List, Set

import pandas as pd
from data_sources import chembl_client


# Allowed short abbreviations
SHORT_WHITELIST: Set[str] = {"bcg", "hcg", "fsh", "lh", "tsh", "gcsf"}


def normalize_name(raw: str) -> List[str]:
    """
    Normalize drug name:
      - Remove parentheses
      - Split combination drugs by '+'
      - Convert lowercase
      - Remove salts (hydrochloride, sulfate, ...)
      - Remove dosage forms (tablet, injection, ...)
      - Keep only letters/numbers/space/hyphen
    """
    if not isinstance(raw, str):
        return []

    name = re.sub(r"\(.*?\)", "", raw)
    parts = re.split(r"\s*\+\s*", name)

    cleaned: List[str] = []

    for p in parts:
        p = p.strip().lower()
        if not p:
            continue

        p = re.sub(
            r"\b(as|as\.)?\s*(hydrochloride|hcl|sulfate|sulphate|na|sodium|"
            r"potassium|phosphate|carbonate|mesilate|mesylate|acetate|tartrate|"
            r"nitrate|maleate|fumarate|succinate|bromide|iodide)\b",
            "",
            p,
        )

        p = re.sub(
            r"\b(tablet|tablets|capsule|capsules|injection|injections|solution|"
            r"suspension|cream|ointment|gel|drops|syrup|patch|spray)\b",
            "",
            p,
        )

        p = re.sub(r"[^a-z0-9\s\-]", "", p)
        p = re.sub(r"\s+", " ", p).strip()

        if p:
            cleaned.append(p)

    return cleaned


def is_radiopharmaceutical(name: str) -> bool:
    """
    Detect radiopharmaceutical patterns:
      - Starting with digits + letters + '-'
      - Patterns like [131i], [32p], ...
    """
    if re.match(r"^\d+\s*[a-z]*-", name):
        return True
    if re.search(r"\[\d+\s*[a-z]+\]", name):
        return True
    return False


def build_normalized_index(names: List[str]) -> Dict[str, List[str]]:
    """
    Build mapping:
      normalized_name -> list of original names
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


def load_local_drugs_from_excel(content: bytes) -> pd.DataFrame:
    """
    Load Excel drug list.
    A column named 'Name' is expected.
    """
    buf = io.BytesIO(content)
    df = pd.read_excel(buf)

    if "Name" not in df.columns:
        raise ValueError("Column 'Name' not found in Excel file.")

    df = df[df["Name"].notna()].copy()
    return df
