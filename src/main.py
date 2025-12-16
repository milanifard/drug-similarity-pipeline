# main.py
from fastapi import FastAPI
from starlette.middleware.sessions import SessionMiddleware
from api.import_local_drugs import router as import_router
from api.similarity import router as similarity_router
from api.sync_chembl import router as chembl_sync_router
from api.check_local_drugs import router as check_router
from src.web import auth
from src.web import dashboard
from src.web.similarity_ui import router as similarity_ui_router
from src.web.import_ui import router as import_ui_router

app = FastAPI(
    title="Drug Similarity API",
    version="2.0"
)

app.add_middleware(
    SessionMiddleware,
    secret_key="OMID_SECRET_MILANIFARD",
    session_cookie="drug_similarity_session",
    same_site="lax",
)

# Register routers
app.include_router(auth.router)
app.include_router(dashboard.router)
app.include_router(similarity_ui_router)
app.include_router(import_ui_router)
app.include_router(check_router)
app.include_router(import_router)
app.include_router(similarity_router)
app.include_router(chembl_sync_router)


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="0.0.0.0", port=8000)
