<div align="center">

# Drug Similarity Pipeline  
### *2D/3D Molecular Similarity Engine using RDKit, ChEMBL, PubChem*

<img src="logo.svg" width="180"/>

---

**Drug Similarity Pipeline (DSP)** is an open-source cheminformatics engine designed to:

✔ Normalize raw drug names from any market  
✔ Map drugs to **ChEMBL** molecules  
✔ Extract **SMILES**  
✔ Build or fetch **3D conformers** (PubChem → RDKit fallback)  
✔ Store everything in a fast local database  
✔ Run **2D Tanimoto** + **3D USRCAT** similarity search  

Built for flexibility — DSP can index any country’s drug list, not only a specific market.

---

</div>

---

# Features

### 🔹 **Name Normalization**
- Removes salt forms (HCl, sulfate, acetate…)
- Removes dosage forms (tablet, injection…)
- Lowercase + cleanup + multi-component splitting
- Radiopharmaceutical filtering
- Works for any country's dataset

---

### 🔹 **ChEMBL Integration**
- Fetches approved-phase drugs  
- Matches normalized names  
- Returns SMILES + ChEMBL IDs  

---

### 🔹 **3D Structure Builder**
Uses a multi-level strategy:

1. **PubChem Bioactive 3D Conformer (preferred)**  
2. **RDKit ETKDGv3 embed** (fallback)

Stores:
- SMILES  
- MolBlock (V2000)  
- Conformer source  

---

### 🔹 **Local Molecule Database**
MariaDB schema includes:

