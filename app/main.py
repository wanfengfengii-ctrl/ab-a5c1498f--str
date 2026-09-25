"""HTTP API 与静态页面入口。"""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse

from .schemas import AdjudicationRequest
from .solver import solve

BASE_DIR = Path(__file__).resolve().parent

app = FastAPI(title="混合 DNA 样本裁决系统", version="1.0.0")


@app.get("/api/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/api/mixture-adjudications")
def adjudicate(request: AdjudicationRequest) -> dict:
    return solve(
        candidates=[
            {"id": c.id, "genotypes": c.genotypes} for c in request.candidates
        ],
        loci=[{"name": l.name, "observed": l.observed} for l in request.loci],
        tolerance=request.tolerance,
    )


@app.get("/", include_in_schema=False)
def index() -> FileResponse:
    return FileResponse(BASE_DIR / "static" / "index.html")
