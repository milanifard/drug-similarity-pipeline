import hashlib
from src.services.chemdb_service import get_engine
from sqlalchemy import text


def hash_password(password: str) -> str:
    return hashlib.sha256(password.encode("utf-8")).hexdigest()


def authenticate_user(username: str, password: str):
    engine = get_engine()
    pwd_hash = hash_password(password)

    with engine.connect() as conn:
        row = conn.execute(
            text("""
                SELECT id, username
                FROM users
                WHERE username = :u
                  AND password_hash = :p
            """),
            {"u": username, "p": pwd_hash}
        ).fetchone()

    return row
