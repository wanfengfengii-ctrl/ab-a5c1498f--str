"""API 层测试：结构校验、状态语义与前端页面可达性。"""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from app.main import app
from tests.cases import (
    ambiguous_payload,
    no_solution_payload,
    unique_payload,
)

client = TestClient(app)


def post(payload: dict):
    return client.post("/api/mixture-adjudications", json=payload)


def test_health():
    resp = client.get("/api/health")
    assert resp.status_code == 200
    assert resp.json() == {"status": "ok"}


def test_index_page_served():
    resp = client.get("/")
    assert resp.status_code == 200
    assert "裁决" in resp.text


def test_unique_via_api():
    resp = post(unique_payload())
    assert resp.status_code == 200
    body = resp.json()
    assert body["status"] == "unique"
    assert body["solution"]["contributors"] == [
        {"id": "C1", "copies": 2},
        {"id": "C3", "copies": 2},
    ]


def test_ambiguous_via_api_lists_two_witnesses():
    resp = post(ambiguous_payload())
    assert resp.status_code == 200
    body = resp.json()
    assert body["status"] == "ambiguous"
    assert body["solution"] is None
    assert len(body["witnesses"]) == 2
    v0 = body["witnesses"][0]["contributors"]
    v1 = body["witnesses"][1]["contributors"]
    assert v0 != v1


def test_no_solution_via_api():
    resp = post(no_solution_payload())
    assert resp.status_code == 200
    assert resp.json()["status"] == "no_solution"


def _too_few_candidates(p):
    p["candidates"] = p["candidates"][:3]


def _too_many_candidates(p):
    for i in range(5, 10):
        p["candidates"].append(
            {"id": f"CX{i}", "genotypes": {l["name"]: ["a", "a"] for l in p["loci"]}}
        )


def _too_few_loci(p):
    p["loci"] = p["loci"][:2]


def _too_many_loci(p):
    for i in range(4, 8):
        p["loci"].append({"name": f"L{i}", "observed": {"a": 4}})


def _duplicate_ids(p):
    p["candidates"][1]["id"] = "C1"


def _inconsistent_depth(p):
    p["loci"][0]["observed"]["a"] = 99


def _negative_tolerance(p):
    p["tolerance"] = -1


def _missing_genotype_locus(p):
    del p["candidates"][0]["genotypes"]["L3"]


def _genotype_not_two_alleles(p):
    p["candidates"][0]["genotypes"]["L1"] = ["a"]


def _negative_observed_count(p):
    p["loci"][0]["observed"]["a"] = -2


def _empty_observed(p):
    p["loci"][0]["observed"] = {}


@pytest.mark.parametrize(
    "mutate",
    [
        _too_few_candidates,
        _too_many_candidates,
        _too_few_loci,
        _too_many_loci,
        _duplicate_ids,
        _inconsistent_depth,
        _negative_tolerance,
        _missing_genotype_locus,
        _genotype_not_two_alleles,
        _negative_observed_count,
        _empty_observed,
    ],
)
def test_invalid_inputs_get_422(mutate):
    payload = unique_payload()
    mutate(payload)
    resp = post(payload)
    assert resp.status_code == 422
    assert resp.json()["detail"]


def test_depth_beyond_limit_gets_422():
    payload = unique_payload()
    for locus in payload["loci"]:
        locus["observed"] = {"a": 501}
    resp = post(payload)
    assert resp.status_code == 422


def test_malformed_json_gets_422():
    resp = client.post(
        "/api/mixture-adjudications",
        content=b"{not-json",
        headers={"Content-Type": "application/json"},
    )
    assert resp.status_code == 422


def test_response_never_disguises_non_unique_as_unique():
    # 歧义与无解响应中不得出现 solution 字段内容
    for payload, status in [(ambiguous_payload(), "ambiguous"), (no_solution_payload(), "no_solution")]:
        body = post(payload).json()
        assert body["status"] == status
        assert body["solution"] is None
