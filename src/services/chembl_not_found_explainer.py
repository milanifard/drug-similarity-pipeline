# src/services/chembl_not_found_explainer.py

import re


def explain_not_in_chembl(normalized_name: str) -> str:
    """
    Return a human-readable explanation why a drug name is not expected
    to be found in ChEMBL.
    """

    name = normalized_name.lower().strip()

    # --------------------------------------------------
    # 1. Radiopharmaceuticals / isotopes
    # --------------------------------------------------
    if re.search(r"\b(tc|technetium|iodine|lutetium|yttrium|gallium|indium|rubidium|samarium|rhenium|radium|thallium)\b", name) \
       or re.search(r"\b\d{2,3}m?\b", name):
        return "Radiopharmaceutical or isotope-based product"

    # --------------------------------------------------
    # 2. Biologics (monoclonal antibodies, proteins, peptides)
    # --------------------------------------------------
    if name.endswith("mab") or name.endswith("umab") or name.endswith("cept"):
        return "Biologic (monoclonal antibody or fusion protein)"

    if re.search(r"\b(insulin|interferon|interleukin|factor\s+(vii|viii|ix|xiii)|epoetin|somatropin|glucagon)\b", name):
        return "Biologic (protein or peptide hormone)"

    # --------------------------------------------------
    # 3. Vaccines / toxoids / antigens
    # --------------------------------------------------
    if re.search(r"\b(vaccine|toxoid|antigen|bcg|measles|mumps|rubella|pertussis|tetanus|diphtheria)\b", name):
        return "Vaccine, toxoid, or antigen-based product"

    # --------------------------------------------------
    # 4. Blood products / cell or tissue therapies
    # --------------------------------------------------
    if re.search(r"\b(albumin|fibrinogen|immunoglobulin|plasma|platelet|keratinocyte|fibroblast|stem\s*cell)\b", name):
        return "Blood product or cell/tissue-based therapy"

    # --------------------------------------------------
    # 5. Inorganic elements, metals, simple salts
    # --------------------------------------------------
    if re.fullmatch(r"(iron|zinc|copper|selenium|sulfur|aluminium|aluminum|calcium|magnesium|potassium|sodium|barium)", name):
        return "Inorganic element or simple mineral salt"

    if re.search(r"\b(chloride|oxide|hydroxide|sulfate|phosphate|carbonate)\b", name):
        return "Inorganic salt or simple chemical compound"

    # --------------------------------------------------
    # 6. Combination drugs / formulations
    # --------------------------------------------------
    if re.search(r"\b(co-|compound|combination|triphasic|biphasic|four-phasic|ld|hd)\b", name):
        return "Combination drug or multi-ingredient formulation"

    if re.search(r"\b(and|plus|with)\b", name):
        return "Combination drug (multiple active ingredients)"

    # --------------------------------------------------
    # 7. Excipients, solvents, non-drug materials
    # --------------------------------------------------
    if re.search(r"\b(ethanol|alcohol|saline|ringer|water|surfactant|glycerin|paraffin|collodion)\b", name):
        return "Excipient, solvent, or non-active pharmaceutical material"

    # --------------------------------------------------
    # 8. Medical procedures / non-molecular concepts
    # --------------------------------------------------
    if re.search(r"\b(burn|dialysis|filtration|haemodialysis|hemofiltration)\b", name):
        return "Medical procedure or non-molecular therapeutic concept"

    # --------------------------------------------------
    # 9. Ambiguous or non-specific names
    # --------------------------------------------------
    if len(name) < 4 or name in {"acid", "drug", "drug name", "digestive"}:
        return "Ambiguous or non-specific drug name"

    # --------------------------------------------------
    # 10. Fallback
    # --------------------------------------------------
    return "Not found in ChEMBL (likely synonym, legacy name, or unsupported category)"
