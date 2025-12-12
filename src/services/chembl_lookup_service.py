from typing import Dict
from src.services.local_drug_service import normalize_name, is_radiopharmaceutical, SHORT_WHITELIST
import pandas as pd


def build_chembl_lookup(df: pd.DataFrame) -> Dict[str, dict]:
    """
    Build lookup table:
        normalized_name -> row from ChEMBL dataframe

    Filtering:
        - length <= 3 (unless whitelisted)
        - skip radiopharmaceuticals
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

            lookup[n] = row

    return lookup
