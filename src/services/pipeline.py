# pipeline.py

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdMolDescriptors

def filter_smiles_for_embedding(df, max_atoms=80):
    allowed = {"H", "C", "N", "O", "F", "S", "Cl", "Br", "I", "P", "B"}
    clean_smiles, clean_names, clean_ids = [], [], []

    for s, n, chembl_id in zip(df["smiles"], df["name"], df["chembl_id"]):
        if not isinstance(s, str) or "." in s:
            continue
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            continue
        if mol.GetNumAtoms() < 3 or mol.GetNumAtoms() > max_atoms:
            continue
        elems = {atom.GetSymbol() for atom in mol.GetAtoms()}
        if not elems.issubset(allowed):
            continue

        clean_smiles.append(s)
        clean_names.append(n)
        clean_ids.append(chembl_id)

    return clean_smiles, clean_names, clean_ids


# ----------------------------
# 3D Similarity (multi confs)
# ----------------------------
def compute_usr_similarity(ref_mol, ref_confs, t_mols, t_confs, names):
    ref_fps = []
    for cid in ref_confs:
        fp = np.array(list(rdMolDescriptors.GetUSRCAT(ref_mol, confId=cid)))
        ref_fps.append(fp)

    results = []
    for mol, conf_ids, name in zip(t_mols, t_confs, names):
        if mol is None or len(conf_ids) == 0:
            results.append((name, 0.0))
            continue

        best = 0.0
        for ref_fp in ref_fps:
            for cid in conf_ids:
                fp = np.array(list(rdMolDescriptors.GetUSRCAT(mol, confId=cid)))
                sim = float(np.dot(ref_fp, fp) /
                            (np.linalg.norm(ref_fp) * np.linalg.norm(fp)))
                if sim > best:
                    best = sim

        results.append((name, best))

    return pd.DataFrame(results, columns=["name", "3D_similarity"])


# ----------------------------
# 2D Similarity
# ----------------------------
def compute_2d_similarity(ref_smiles, smiles_list, names, radius=2, nbits=2048):
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, radius, nBits=nbits)

    rows = []
    for s, name in zip(smiles_list, names):
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            rows.append((name, 0.0))
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
        sim = DataStructs.TanimotoSimilarity(ref_fp, fp)
        rows.append((name, sim))

    return pd.DataFrame(rows, columns=["name", "2D_similarity"])
