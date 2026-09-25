"""求解器单元测试：直接调用 solve()，不经过 HTTP 层。"""

from __future__ import annotations

from app.solver import solve
from tests.cases import (
    ambiguous_payload,
    fewer_contributors_win_payload,
    no_solution_payload,
    unique_payload,
)


def run(payload: dict) -> dict:
    return solve(
        candidates=payload["candidates"],
        loci=payload["loci"],
        tolerance=payload["tolerance"],
    )


def test_unique_solution_recovers_true_mixture():
    res = run(unique_payload())
    assert res["status"] == "unique"
    assert res["objective"] == {"contributors": 2, "absoluteResidual": 0}
    assert res["solution"]["contributors"] == [
        {"id": "C1", "copies": 2},
        {"id": "C3", "copies": 2},
    ]
    assert res["totalDepth"] == 4
    # 可复核：逐位点预测与观测完全一致
    for locus in res["solution"]["perLocus"]:
        for row in locus["alleles"]:
            assert row["predicted"] == row["observed"]
            assert row["absDiff"] == 0


def test_ambiguous_when_two_vectors_share_optimum():
    res = run(ambiguous_payload())
    assert res["status"] == "ambiguous"
    assert res["solution"] is None
    witnesses = res["witnesses"]
    assert len(witnesses) == 2
    vectors = [tuple((c["id"], c["copies"]) for c in w["contributors"]) for w in witnesses]
    assert vectors[0] != vectors[1]
    assert vectors == [(("C1", 4),), (("C2", 4),)]
    assert res["objective"] == {"contributors": 1, "absoluteResidual": 0}


def test_fewer_contributors_beat_lower_residual():
    res = run(fewer_contributors_win_payload())
    assert res["status"] == "unique"
    # k=1 残差 6（每位点 |4-3|+|0-1|=2，共 3 个位点）优于任何 k=2 的零残差组合（人数优先）
    assert res["objective"] == {"contributors": 1, "absoluteResidual": 6}
    assert res["solution"]["contributors"] == [{"id": "C1", "copies": 4}]


def test_no_solution_is_reported_plainly():
    res = run(no_solution_payload())
    assert res["status"] == "no_solution"
    assert res["objective"] is None
    assert res["solution"] is None
    assert res["witnesses"] is None


def test_unobserved_alleles_are_penalized():
    # C1 携带观测中不存在的等位基因 z，残差必须计入 |0 - 预测|
    payload = {
        "tolerance": 2,
        "candidates": [
            {"id": "C1", "genotypes": {n: ["a", "z"] for n in ["L1", "L2", "L3"]}},
            {"id": "C2", "genotypes": {n: ["b", "b"] for n in ["L1", "L2", "L3"]}},
            {"id": "C3", "genotypes": {n: ["c", "c"] for n in ["L1", "L2", "L3"]}},
            {"id": "C4", "genotypes": {n: ["d", "d"] for n in ["L1", "L2", "L3"]}},
        ],
        "loci": [
            {"name": n, "observed": {"a": 2}} for n in ["L1", "L2", "L3"]
        ],
    }
    res = run(payload)
    assert res["status"] == "unique"
    # C1×2：每位点 pred a=1、z=1 → 残差 |1-2|+|1-0|=2，三个位点共 6
    assert res["objective"] == {"contributors": 1, "absoluteResidual": 6}
    locus = res["solution"]["perLocus"][0]
    zrow = next(r for r in locus["alleles"] if r["allele"] == "z")
    assert zrow["observed"] == 0 and zrow["predicted"] == 1 and zrow["absDiff"] == 1


def test_heterozygote_contributes_half_copies():
    # 单人杂合 4 份 → 每条等位基因预测 2 份（半份分摊模型）
    payload = {
        "tolerance": 0,
        "candidates": [
            {"id": "C1", "genotypes": {n: ["a", "b"] for n in ["L1", "L2", "L3"]}},
            {"id": "C2", "genotypes": {n: ["c", "c"] for n in ["L1", "L2", "L3"]}},
            {"id": "C3", "genotypes": {n: ["d", "d"] for n in ["L1", "L2", "L3"]}},
            {"id": "C4", "genotypes": {n: ["e", "e"] for n in ["L1", "L2", "L3"]}},
        ],
        "loci": [{"name": n, "observed": {"a": 2, "b": 2}} for n in ["L1", "L2", "L3"]],
    }
    res = run(payload)
    assert res["status"] == "unique"
    assert res["solution"]["contributors"] == [{"id": "C1", "copies": 4}]


def test_three_contributor_solution_found():
    # 构造必须 3 人才能解释的样本：C1×1 + C2×1 + C3×2，总深度 4
    loci_names = ["L1", "L2", "L3"]
    payload = {
        "tolerance": 0,
        "candidates": [
            {"id": "C1", "genotypes": {n: ["a", "a"] for n in loci_names}},
            {"id": "C2", "genotypes": {n: ["b", "b"] for n in loci_names}},
            {"id": "C3", "genotypes": {n: ["c", "d"] for n in loci_names}},
            {"id": "C4", "genotypes": {n: ["e", "e"] for n in loci_names}},
        ],
        "loci": [
            {"name": n, "observed": {"a": 1, "b": 1, "c": 1, "d": 1}}
            for n in loci_names
        ],
    }
    res = run(payload)
    assert res["status"] == "unique"
    assert res["objective"] == {"contributors": 3, "absoluteResidual": 0}
    assert res["solution"]["contributors"] == [
        {"id": "C1", "copies": 1},
        {"id": "C2", "copies": 1},
        {"id": "C3", "copies": 2},
    ]


def test_examined_counters_are_reported():
    res = run(unique_payload())
    assert res["examined"]["subsets"] > 0
    assert res["examined"]["assignments"] > 0
