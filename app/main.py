"""FastAPI 入口：静态页面 + POST /api/mixture-adjudications。"""

from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from .models import AdjudicationRequest
from .solver import solve

app = FastAPI(title="混合 DNA 样本供样人裁决服务")

STATIC_DIR = Path(__file__).parent / "static"


@app.get("/api/health")
def health():
    return {"status": "ok"}


@app.post("/api/mixture-adjudications")
def adjudicate(request: AdjudicationRequest):
    genotypes = {c.id: c.genotypes for c in request.candidates}
    return solve(
        [c.id for c in request.candidates],
        genotypes,
        request.loci,
        request.observed,
        request.tolerance,
    )


@app.get("/", include_in_schema=False)
def index():
    return FileResponse(STATIC_DIR / "index.html")


app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
