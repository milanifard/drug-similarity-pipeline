<div align="center">

# Drug Similarity Pipeline  
### *2D/3D Molecular Similarity Engine using RDKit, ChEMBL, PubChem*

<img src="logo.svg" width="180"/>

---

Drug Similarity Pipeline

Drug Similarity Pipeline is an open-source system for:

Normalizing drug names from national market lists (Excel files)

Mapping normalized names to approved molecules in ChEMBL

Extracting canonical SMILES

Building or fetching 3D conformers (PubChem → RDKit fallback)

Storing all results in a MariaDB database

Performing similarity search using combined 2D (Tanimoto) and 3D (USRCAT) scores

The entire project is fully Dockerized and runs through docker-compose.

Features

Drug name normalization:

Remove salts, dosage forms, parentheses, and noise tokens

Split combination drugs (A + B)

Lowercase normalization

Remove radiopharmaceuticals

Matching normalized drug names to approved ChEMBL molecules

Extracting SMILES and building 3D conformers:

Fetch PubChem bioactive conformer if available

Otherwise generate ETKDG 3D conformers using RDKit

Storing results in a local MariaDB table (drug_molecules or your preferred name)

normalized_name

original_name

chembl_name

chembl_id

smiles

molblock_3d (V2000 MolBlock)

conformer_source (pubchem or rdkit)

Similarity search API endpoint using only the locally stored drug dataset

High-Level Architecture

Components:

FastAPI service (image: drug-similarity-api)

Normalization, ChEMBL lookup, PubChem retrieval

Conformer generation

Similarity computation (2D + 3D)

MariaDB service (image: chemdb)

Stores normalized drug records and conformers

Communication occurs entirely through Docker’s internal network.
Database access string is passed via the environment variable CHEMDB_URL.

Requirements

Docker

docker-compose

Internet access (required for ChEMBL and PubChem)

Running with Docker

Clone the repository:
git clone https://github.com/milanifard/drug-similarity-pipeline.git
cd drug-similarity-pipeline


    Check docker-compose.yml.

It should define:

    drug-similarity-api service

    chemdb MariaDB service

    Environment variable such as:

CHEMDB_URL=mysql+pymysql://chemuser:chempass@chemdb:3306/chemdb

    Start the system:

docker-compose up --build

This will:

    Start MariaDB

    Start the FastAPI service

    Automatically run any SQL initialization scripts (if mounted under /docker-entrypoint-initdb.d)

    API URL:

http://localhost:8000

Swagger UI:

http://localhost:8000/docs

Database Initialization

Place schema files such as:

db/init/01_create_schema.sql

and mount them in docker-compose.yml under:

/docker-entrypoint-initdb.d

This ensures the required tables are created automatically on first run.

To reset the database fully:

docker-compose down -v
docker-compose up --build

API Endpoints
1. Import and Normalize Drug List

POST /import_drugs

Input:

    Excel file (containing a column of drug names)

Internal workflow:
    Load all drug names from the Excel file
    Normalize names
    Fetch approved ChEMBL molecules
    Match normalized names to ChEMBL

    For each new drug:
        Retrieve SMILES
        Fetch or generate 3D conformer
        Store into database (drug_molecules)
    Log progress in the API container console (processed, skipped, existing, etc.)

Example output:

{
  "inserted": 350,
  "exists": 120,
  "skipped": 800
}

2. Similarity Search

GET /similarity

Parameters:

    drug: reference drug name

    alpha: weight for combining 2D vs 3D similarity scores

Workflow:

    Look up reference drug in ChEMBL

    Retrieve or generate reference conformer

    Load all target molecules from local database (drug_molecules)

    Compute:

        2D Tanimoto similarity

        3D USRCAT similarity

    Produce a ranked similarity list

Example output:

{
  "query": "metformin",
  "results": [
    {
      "name": "METFORMIN HYDROCHLORIDE",
      "normalized_name": "metformin",
      "chembl_id": "CHEMBL1431",
      "similarity_2d": 0.87,
      "similarity_3d": 0.79,
      "weighted_similarity": 0.82
    }
  ]
}

Project Structure

src/
  api.py
  pipeline.py
  chemdb_service.py
  drug_import_service.py
  pubchem_utils.py
  data_sources/
    chembl_client.py

docker-compose.yml
Dockerfile
requirements.txt

License

MIT License (or any license you prefer).
Add a LICENSE file at the root of the project