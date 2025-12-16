from fastapi import APIRouter, Request, Form
from fastapi.responses import RedirectResponse
from fastapi.templating import Jinja2Templates
from src.services.auth_service import authenticate_user

router = APIRouter()
templates = Jinja2Templates(directory="src/templates")


@router.get("/login")
def login_page(request: Request):
    return templates.TemplateResponse(
        "login.html",
        {"request": request}
    )


@router.post("/login")
def login_action(
    request: Request,
    username: str = Form(...),
    password: str = Form(...)
):
    user = authenticate_user(username, password)
    if not user:
        return RedirectResponse("/login?error=1", status_code=302)

    request.session["user"] = {
        "id": user.id,
        "username": user.username,
    }

    return RedirectResponse("/dashboard", status_code=302)

@router.get("/logout")
def logout(request: Request):
    request.session.clear()
    return RedirectResponse(url="/login", status_code=302)
