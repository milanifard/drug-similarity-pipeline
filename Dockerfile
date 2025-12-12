# ============================================================
# Dockerfile for FastAPI + RDKit molecular similarity API
# Base: micromamba (Jammy), RDKit via conda-forge
# ============================================================

FROM mambaorg/micromamba:1.5.8-jammy

# همه‌چیز در env پایه نصب می‌شود
ARG MAMBA_DOCKERFILE_ACTIVATE=1
WORKDIR /app

# RDKit + numpy + pandas + requests از conda-forge
# (پین نسخه‌ها برای سازگاری بهتر)
RUN micromamba install -y -n base -c conda-forge \
      python=3.10 \
      rdkit=2023.09.5 \
      numpy=1.26 \
      pandas=2.2 \
      requests \
    && micromamba clean -a -y

# بسته‌های pip (FastAPI و …)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# کدها
COPY ./src /app/src

# پورت FastAPI
EXPOSE 8000
ENV PYTHONPATH=/app/src
# اجرای سرویس
CMD ["uvicorn", "src.main:app", "--host", "0.0.0.0", "--port", "8000"]
