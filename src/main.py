# main.py
from fastapi import FastAPI
from api.import_local_drugs import router as import_router
from api.similarity import router as similarity_router
from api.sync_chembl import router as chembl_sync_router

app = FastAPI(
    title="Drug Similarity API",
    version="2.0"
)

# Register routers
app.include_router(import_router)
app.include_router(similarity_router)
app.include_router(chembl_sync_router)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="0.0.0.0", port=8000)
