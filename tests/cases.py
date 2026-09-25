"""共享的测试用例构造。"""

from __future__ import annotations


def unique_payload() -> dict:
    """真实混合物为 C1×2 + C3×2（总深度 4），误差 0 下应得唯一解。"""
    return {
        "tolerance": 0,
        "candidates": [
            {"id": "C1", "genotypes": {"L1": ["a", "b"], "L2": ["x", "y"], "L3": ["p", "q"]}},
            {"id": "C2", "genotypes": {"L1": ["a", "c"], "L2": ["x", "x"], "L3": ["p", "p"]}},
            {"id": "C3", "genotypes": {"L1": ["b", "c"], "L2": ["y", "z"], "L3": ["q", "r"]}},
            {"id": "C4", "genotypes": {"L1": ["c", "c"], "L2": ["z", "z"], "L3": ["r", "r"]}},
        ],
        "loci": [
            {"name": "L1", "observed": {"a": 1, "b": 2, "c": 1}},
            {"name": "L2", "observed": {"x": 1, "y": 2, "z": 1}},
            {"name": "L3", "observed": {"p": 1, "q": 2, "r": 1}},
        ],
    }


def ambiguous_payload() -> dict:
    """C1 与 C2 基因型完全相同，单人 4 份即可精确解释观测 → 歧义。"""
    geno = {"L1": ["a", "b"], "L2": ["x", "y"], "L3": ["p", "q"]}
    return {
        "tolerance": 0,
        "candidates": [
            {"id": "C1", "genotypes": dict(geno)},
            {"id": "C2", "genotypes": dict(geno)},
            {"id": "C3", "genotypes": {"L1": ["a", "c"], "L2": ["x", "z"], "L3": ["p", "r"]}},
            {"id": "C4", "genotypes": {"L1": ["c", "c"], "L2": ["z", "z"], "L3": ["r", "r"]}},
        ],
        "loci": [
            {"name": "L1", "observed": {"a": 2, "b": 2}},
            {"name": "L2", "observed": {"x": 2, "y": 2}},
            {"name": "L3", "observed": {"p": 2, "q": 2}},
        ],
    }


def fewer_contributors_win_payload() -> dict:
    """k=1 可行（残差 6），尽管 k=2 存在零残差组合，也必须先最少化人数。"""
    loci_names = ["L1", "L2", "L3"]
    return {
        "tolerance": 1,
        "candidates": [
            {"id": "C1", "genotypes": {n: ["a", "a"] for n in loci_names}},
            {"id": "C2", "genotypes": {"L1": ["a", "b"], "L2": ["c", "c"], "L3": ["c", "c"]}},
            {"id": "C3", "genotypes": {n: ["a", "c"] for n in loci_names}},
            {"id": "C4", "genotypes": {n: ["d", "d"] for n in loci_names}},
        ],
        "loci": [{"name": n, "observed": {"a": 3, "b": 1}} for n in loci_names],
    }


def no_solution_payload() -> dict:
    """观测中出现任何候选人都不具备的等位基因 zz，误差 0 下无解。"""
    payload = unique_payload()
    payload["loci"][0]["observed"] = {"a": 1, "b": 2, "zz": 1}
    return payload
