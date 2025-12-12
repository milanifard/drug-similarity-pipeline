# src/services/chembl_lookup_service.py

from typing import Dict
from src.services.drug_name_normalizer import normalize_name, is_radiopharmaceutical, SHORT_WHITELIST, is_banned_term
import pandas as pd


def build_chembl_lookup(df: pd.DataFrame) -> Dict[str, dict]:
    """
    Build lookup mapping:
        normalized_name -> row from ChEMBL dataframe

    Conservative normalization:
        - skip radiopharmaceuticals
        - skip very short names unless whitelisted
        - skip banned non-drug categories
    """

    lookup = {}

    for _, row in df.iterrows():
        raw_name = row["name"]
        norm_list = normalize_name(raw_name)

        if not norm_list:
            continue

        for n in norm_list:
            if len(n) <= 3 and n not in SHORT_WHITELIST:
                continue

            if is_radiopharmaceutical(n):
                continue

            if is_banned_term(n):
                continue

            lookup[n] = row

    return lookup
