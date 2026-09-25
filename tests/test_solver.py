"""求解器单元测试：唯一解、歧义、无解、目标优先级、纯合子、半整数残差。"""

from app.solver import solve

LOCI = ["L1", "L2", "L3"]

# 内置示例：真实混合物为 C1×4 + C3×2，总深度 6，误差上限 0 时有唯一解。
EXAMPLE_GENOTYPES = {
    "C1": {"L1": ("12", "13"), "L2": ("7", "9"), "L3": ("16", "17")},
    "C2": {"L1": ("11", "12"), "L2": ("8", "8"), "L3": ("15", "18")},
    "C3": {"L1": ("14", "14"), "L2": ("9", "10"), "L3": ("17", "18")},
    "C4": {"L1": ("13", "15"), "L2": ("7", "10"), "L3": ("16", "16")},
    "C5": {"L1": ("11", "15"), "L2": ("8", "11"), "L3": ("15", "15")},
}
EXAMPLE_OBSERVED = {
    "L1": {"12": 2, "13": 2, "14": 2},
    "L2": {"7": 2, "9": 3, "10": 1},
    "L3": {"16": 2, "17": 3, "18": 1},
}


def test_unique_example():
    r = solve(list(EXAMPLE_GENOTYPES), EXAMPLE_GENOTYPES, LOCI, EXAMPLE_OBSERVED, 0)
    assert r["status"] == "unique"
    assert r["totalDepth"] == 6
    assert r["objective"] == {"contributors": 2, "residualSum": 0}
    assert r["optimalVectorCount"] == 1
    assert r["solution"]["contributors"] == [
        {"candidateId": "C1", "copies": 4},
        {"candidateId": "C3", "copies": 2},
    ]
    # 逐位点复核：所有残差为 0
    for locus in r["solution"]["loci"]:
        for row in locus["alleles"]:
            assert row["residual"] == 0
            assert row["predicted"] == row["observed"]


def test_unique_with_noise_within_tolerance():
    # L1 观测噪声：12 号多 1 份、13 号少 1 份（总深度仍为 6）
    # 最优解不变，残差和为 2，且仍在误差上限 1 内
    observed = {l: dict(a) for l, a in EXAMPLE_OBSERVED.items()}
    observed["L1"]["12"] = 3
    observed["L1"]["13"] = 1
    r = solve(list(EXAMPLE_GENOTYPES), EXAMPLE_GENOTYPES, LOCI, observed, 1)
    assert r["status"] == "unique"
    assert r["objective"]["contributors"] == 2
    assert r["objective"]["residualSum"] == 2
    assert r["solution"]["contributors"] == [
        {"candidateId": "C1", "copies": 4},
        {"candidateId": "C3", "copies": 2},
    ]


def test_ambiguous_identical_genotypes():
    # P1 与 P2 基因型完全相同：P1×4 与 P2×4 同为最优 → 歧义
    genotypes = {
        "P1": {"L1": ("9", "10"), "L2": ("11", "12"), "L3": ("13", "14")},
        "P2": {"L1": ("9", "10"), "L2": ("11", "12"), "L3": ("13", "14")},
        "P3": {"L1": ("8", "8"), "L2": ("8", "8"), "L3": ("8", "8")},
        "P4": {"L1": ("7", "7"), "L2": ("7", "7"), "L3": ("7", "7")},
    }
    observed = {
        "L1": {"9": 2, "10": 2},
        "L2": {"11": 2, "12": 2},
        "L3": {"13": 2, "14": 2},
    }
    r = solve(list(genotypes), genotypes, LOCI, observed, 0)
    assert r["status"] == "ambiguous"
    assert r["objective"] == {"contributors": 1, "residualSum": 0}
    assert r["optimalVectorCount"] == 2
    assert r["solution"] is None
    assert len(r["witnesses"]) == 2
    w0, w1 = r["witnesses"]
    assert w0["contributors"] != w1["contributors"]
    assert {w0["contributors"][0]["candidateId"], w1["contributors"][0]["candidateId"]} == {"P1", "P2"}
    assert w0["residualSum"] == w1["residualSum"] == 0


def test_ambiguous_two_pairs():
    # 两对不同供样人都能完美解释样本 → 歧义，见证分属两对
    genotypes = {
        "P1": {"L1": ("9", "10"), "L2": ("13", "14"), "L3": ("17", "18")},
        "P2": {"L1": ("11", "12"), "L2": ("15", "16"), "L3": ("19", "20")},
        "P3": {"L1": ("9", "11"), "L2": ("13", "15"), "L3": ("17", "19")},
        "P4": {"L1": ("10", "12"), "L2": ("14", "16"), "L3": ("18", "20")},
    }
    observed = {
        "L1": {"9": 1, "10": 1, "11": 1, "12": 1},
        "L2": {"13": 1, "14": 1, "15": 1, "16": 1},
        "L3": {"17": 1, "18": 1, "19": 1, "20": 1},
    }
    r = solve(list(genotypes), genotypes, LOCI, observed, 0)
    assert r["status"] == "ambiguous"
    assert r["objective"] == {"contributors": 2, "residualSum": 0}
    assert r["optimalVectorCount"] == 2
    pairs = [tuple(c["candidateId"] for c in w["contributors"]) for w in r["witnesses"]]
    assert sorted(pairs) == [("P1", "P2"), ("P3", "P4")]
    assert all(c["copies"] == 2 for w in r["witnesses"] for c in w["contributors"])


def test_no_solution():
    # 观测等位基因 99 不在任何候选人基因型中，误差上限 0 → 无解
    observed = {l: dict(a) for l, a in EXAMPLE_OBSERVED.items()}
    observed["L1"] = {"99": 6}
    r = solve(list(EXAMPLE_GENOTYPES), EXAMPLE_GENOTYPES, LOCI, observed, 0)
    assert r["status"] == "no_solution"
    assert r["objective"] is None
    assert r["solution"] is None
    assert r["witnesses"] == []


def test_contributor_count_minimized_before_residual():
    # 单人解释（残差和 1）与双人完美解释（残差和 0）同时可行时，
    # 必须先最少化供样人数 → 选单人方案
    genotypes = {
        "P1": {"L1": ("9", "9"), "L2": ("9", "9"), "L3": ("9", "9")},
        "P2": {"L1": ("10", "10"), "L2": ("10", "10"), "L3": ("10", "10")},
        "P3": {"L1": ("9", "10"), "L2": ("9", "10"), "L3": ("9", "10")},
        "P4": {"L1": ("7", "7"), "L2": ("7", "7"), "L3": ("7", "7")},
    }
    observed = {
        "L1": {"9": 2, "10": 3},
        "L2": {"9": 2, "10": 3},
        "L3": {"9": 2, "10": 3},
    }
    # P3×5 预测 9:2.5 / 10:2.5（每位点残差和 1，三个位点共 3，单人可行）；
    # P1×2 + P2×3 预测 9:2 / 10:3（残差和 0，双人可行）。
    r = solve(list(genotypes), genotypes, LOCI, observed, 1)
    assert r["status"] == "unique"
    assert r["objective"] == {"contributors": 1, "residualSum": 3}
    assert r["solution"]["contributors"] == [{"candidateId": "P3", "copies": 5}]


def test_homozygote_dosage():
    # 纯合子剂量为 2：P1×3 在每位点贡献 9 号等位基因 3 份
    genotypes = {
        "P1": {"L1": ("9", "9"), "L2": ("9", "9"), "L3": ("9", "9")},
        "P2": {"L1": ("10", "11"), "L2": ("10", "11"), "L3": ("10", "11")},
        "P3": {"L1": ("12", "13"), "L2": ("12", "13"), "L3": ("12", "13")},
        "P4": {"L1": ("14", "15"), "L2": ("14", "15"), "L3": ("14", "15")},
    }
    observed = {"L1": {"9": 3}, "L2": {"9": 3}, "L3": {"9": 3}}
    r = solve(list(genotypes), genotypes, LOCI, observed, 0)
    assert r["status"] == "unique"
    assert r["solution"]["contributors"] == [{"candidateId": "P1", "copies": 3}]
    for locus in r["solution"]["loci"]:
        assert locus["alleles"] == [
            {"allele": "9", "observed": 3, "predicted": 3, "residual": 0}
        ]


def test_half_integer_prediction_and_residual():
    # 奇数份数的杂合子产生 x.5 预测；残差和按真实值累计
    genotypes = {
        "P1": {"L1": ("9", "10"), "L2": ("9", "10"), "L3": ("9", "10")},
        "P2": {"L1": ("9", "9"), "L2": ("9", "9"), "L3": ("9", "9")},
        "P3": {"L1": ("10", "10"), "L2": ("10", "10"), "L3": ("10", "10")},
        "P4": {"L1": ("9", "11"), "L2": ("9", "11"), "L3": ("9", "11")},
    }
    observed = {
        "L1": {"9": 2, "10": 1},
        "L2": {"9": 2, "10": 1},
        "L3": {"9": 2, "10": 1},
    }
    # P1×3：预测 9:1.5 / 10:1.5 → 每位点残差 0.5+0.5=1，三个位点共 3
    r = solve(list(genotypes), genotypes, LOCI, observed, 1)
    assert r["status"] == "unique"
    assert r["objective"] == {"contributors": 1, "residualSum": 3}
    locus = r["solution"]["loci"][0]
    rows = {a["allele"]: a for a in locus["alleles"]}
    assert rows["9"]["predicted"] == 1.5 and rows["9"]["residual"] == -0.5
    assert rows["10"]["predicted"] == 1.5 and rows["10"]["residual"] == 0.5


def test_tolerance_boundary():
    # 残差恰等于误差上限时可行；上限不足则无解
    genotypes = {
        "P1": {"L1": ("9", "9"), "L2": ("9", "9"), "L3": ("9", "9")},
        "P2": {"L1": ("10", "11"), "L2": ("10", "11"), "L3": ("10", "11")},
        "P3": {"L1": ("12", "13"), "L2": ("12", "13"), "L3": ("12", "13")},
        "P4": {"L1": ("14", "15"), "L2": ("14", "15"), "L3": ("14", "15")},
    }
    observed = {"L1": {"9": 4}, "L2": {"9": 4}, "L3": {"9": 4}}
    r = solve(list(genotypes), genotypes, LOCI, observed, 0)
    assert r["status"] == "unique"
    assert r["solution"]["contributors"] == [{"candidateId": "P1", "copies": 4}]

    # L1 观测混入 1 份无人携带的等位基因 10（9 相应少 1 份，总深度仍为 4）：
    # P1×4 的残差为 9 号 1 份 + 10 号 1 份 = 2，恰在上限 1 内逐等位基因可行
    observed_off = {"L1": {"9": 3, "10": 1}, "L2": {"9": 4}, "L3": {"9": 4}}
    r1 = solve(list(genotypes), genotypes, LOCI, observed_off, 1)
    assert r1["status"] == "unique"
    assert r1["objective"] == {"contributors": 1, "residualSum": 2}
    # 上限 0 时同样的观测不可行 → 无解
    r0 = solve(list(genotypes), genotypes, LOCI, observed_off, 0)
    assert r0["status"] == "no_solution"


def test_depth_one_single_contributor():
    genotypes = {
        "P1": {"L1": ("9", "9"), "L2": ("9", "9"), "L3": ("9", "9")},
        "P2": {"L1": ("10", "11"), "L2": ("10", "11"), "L3": ("10", "11")},
        "P3": {"L1": ("12", "13"), "L2": ("12", "13"), "L3": ("12", "13")},
        "P4": {"L1": ("14", "15"), "L2": ("14", "15"), "L3": ("14", "15")},
    }
    observed = {"L1": {"9": 1}, "L2": {"9": 1}, "L3": {"9": 1}}
    r = solve(list(genotypes), genotypes, LOCI, observed, 0)
    assert r["status"] == "unique"
    assert r["solution"]["contributors"] == [{"candidateId": "P1", "copies": 1}]
