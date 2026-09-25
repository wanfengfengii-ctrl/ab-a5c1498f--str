"""API 集成测试：真实 HTTP 往返（FastAPI TestClient）。"""

import copy

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

VALID_PAYLOAD = {
    "candidates": [
        {"id": "C1", "genotypes": {"L1": ["12", "13"], "L2": ["7", "9"], "L3": ["16", "17"]}},
        {"id": "C2", "genotypes": {"L1": ["11", "12"], "L2": ["8", "8"], "L3": ["15", "18"]}},
        {"id": "C3", "genotypes": {"L1": ["14", "14"], "L2": ["9", "10"], "L3": ["17", "18"]}},
        {"id": "C4", "genotypes": {"L1": ["13", "15"], "L2": ["7", "10"], "L3": ["16", "16"]}},
        {"id": "C5", "genotypes": {"L1": ["11", "15"], "L2": ["8", "11"], "L3": ["15", "15"]}},
    ],
    "loci": ["L1", "L2", "L3"],
    "observed": {
        "L1": {"12": 2, "13": 2, "14": 2},
        "L2": {"7": 2, "9": 3, "10": 1},
        "L3": {"16": 2, "17": 3, "18": 1},
    },
    "tolerance": 0,
}


def post(payload):
    return client.post("/api/mixture-adjudications", json=payload)


def test_health():
    r = client.get("/api/health")
    assert r.status_code == 200
    assert r.json() == {"status": "ok"}


def test_index_page_served():
    r = client.get("/")
    assert r.status_code == 200
    assert "裁决" in r.text


def test_static_assets_served():
    for path in ("/static/app.js", "/static/style.css"):
        r = client.get(path)
        assert r.status_code == 200
        assert len(r.content) > 100


def test_unique_adjudication():
    r = post(VALID_PAYLOAD)
    assert r.status_code == 200
    body = r.json()
    assert body["status"] == "unique"
    assert body["totalDepth"] == 6
    assert body["objective"] == {"contributors": 2, "residualSum": 0}
    assert body["solution"]["contributors"] == [
        {"candidateId": "C1", "copies": 4},
        {"candidateId": "C3", "copies": 2},
    ]
    assert body["witnesses"] is None
    assert len(body["solution"]["loci"]) == 3


def test_ambiguous_adjudication():
    payload = {
        "candidates": [
            {"id": "P1", "genotypes": {"L1": ["9", "10"], "L2": ["11", "12"], "L3": ["13", "14"]}},
            {"id": "P2", "genotypes": {"L1": ["9", "10"], "L2": ["11", "12"], "L3": ["13", "14"]}},
            {"id": "P3", "genotypes": {"L1": ["8", "8"], "L2": ["8", "8"], "L3": ["8", "8"]}},
            {"id": "P4", "genotypes": {"L1": ["7", "7"], "L2": ["7", "7"], "L3": ["7", "7"]}},
        ],
        "loci": ["L1", "L2", "L3"],
        "observed": {
            "L1": {"9": 2, "10": 2},
            "L2": {"11": 2, "12": 2},
            "L3": {"13": 2, "14": 2},
        },
        "tolerance": 0,
    }
    r = post(payload)
    assert r.status_code == 200
    body = r.json()
    assert body["status"] == "ambiguous"
    assert body["solution"] is None
    assert body["optimalVectorCount"] == 2
    assert len(body["witnesses"]) == 2
    assert body["witnesses"][0]["contributors"] != body["witnesses"][1]["contributors"]


def test_no_solution_adjudication():
    payload = copy.deepcopy(VALID_PAYLOAD)
    payload["observed"]["L1"] = {"99": 6}
    r = post(payload)
    assert r.status_code == 200
    body = r.json()
    assert body["status"] == "no_solution"
    assert body["solution"] is None
    assert body["witnesses"] == []


def test_unequal_depths_rejected():
    payload = copy.deepcopy(VALID_PAYLOAD)
    payload["observed"]["L2"]["7"] = 5  # L2 总深度变为 9 ≠ 6
    r = post(payload)
    assert r.status_code == 422
    assert "总深度" in str(r.json())


def test_too_few_candidates_rejected():
    payload = copy.deepcopy(VALID_PAYLOAD)
    payload["candidates"] = payload["candidates"][:3]
    r = post(payload)
    assert r.status_code == 422


def test_too_many_candidates_rejected():
    payload = copy.deepcopy(VALID_PAYLOAD)
    for i in range(6, 11):
        payload["candidates"].append(
            {"id": f"C{i}", "genotypes": {"L1": ["1", "1"], "L2": ["1", "1"], "L3": ["1", "1"]}}
        )
    r = post(payload)
    assert r.status_code == 422


def test_duplicate_candidate_ids_rejected():
    payload = copy.deepcopy(VALID_PAYLOAD)
    payload["candidates"][1]["id"] = "C1"
    r = post(payload)
    assert r.status_code == 422
    assert "唯一" in str(r.json())


def test_missing_genotype_locus_rejected():
    payload = copy.deepcopy(VALID_PAYLOAD)
    del payload["candidates"][0]["genotypes"]["L3"]
    r = post(payload)
    assert r.status_code == 422


def test_wrong_allele_count_rejected():
    payload = copy.deepcopy(VALID_PAYLOAD)
    payload["candidates"][0]["genotypes"]["L1"] = ["12"]
    r = post(payload)
    assert r.status_code == 422


def test_negative_tolerance_rejected():
    payload = copy.deepcopy(VALID_PAYLOAD)
    payload["tolerance"] = -1
    r = post(payload)
    assert r.status_code == 422


def test_negative_observed_count_rejected():
    payload = copy.deepcopy(VALID_PAYLOAD)
    payload["observed"]["L1"]["12"] = -2
    r = post(payload)
    assert r.status_code == 422


def test_missing_observed_locus_rejected():
    payload = copy.deepcopy(VALID_PAYLOAD)
    del payload["observed"]["L2"]
    r = post(payload)
    assert r.status_code == 422


def test_too_few_loci_rejected():
    payload = copy.deepcopy(VALID_PAYLOAD)
    payload["loci"] = ["L1", "L2"]
    for c in payload["candidates"]:
        del c["genotypes"]["L3"]
    del payload["observed"]["L3"]
    r = post(payload)
    assert r.status_code == 422
