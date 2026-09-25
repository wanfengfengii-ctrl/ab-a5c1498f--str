"""混合 STR 样本供样人裁决核心求解器。

模型定义
--------
每位候选供样人在每个位点具有两条等位基因（可以相同，即纯合子）。
若供样人 i 被选中且模板份数为 c_i（正整数），则其在某位点对等位基因 a 的
预测拷贝数贡献为 c_i × (a 在其基因型中的剂量 0/1/2) / 2。
因此每个位点的预测总深度恰等于所有选中者份数之和，而该和必须等于各位点
共同的观测总深度 D。

约束：对每个位点、每种相关等位基因（观测中出现、或选中者基因型中出现），
|预测拷贝数 - 观测拷贝数| ≤ 统一误差上限。

目标：先最少化供样人数（1~3 人），再最小化全部绝对残差之和。

唯一性：达到同一最优值的“按候选人编号排序的份数向量”若只有一种，则为唯一
解；否则判定为歧义并给出两份不同见证。无解与歧义都不会被伪装成唯一解。

内部实现：为避免半整数浮点误差，所有预测与观测都乘以 2 后在整数域比较
（|2·预测 - 2·观测| ≤ 2·误差上限），展示时再除回。
"""

from __future__ import annotations

from itertools import combinations
from typing import Dict, List, Sequence, Tuple

import numpy as np

MAX_DEPTH = 1000           # 单个位点允许的观测总深度上限
MAX_OBSERVED_ALLELES = 32  # 单个位点允许的观测等位基因种数上限
_BATCH = 50_000            # 份数组合的分块大小

# 份数向量：按候选人编号排序的 ((编号, 份数), ...) 元组
Vector = Tuple[Tuple[str, int], ...]


def _natural_key(text: str):
    """数字编号按数值排序，其余按字典序，保证确定性。"""
    s = text.strip()
    return (0, int(s)) if s.isdigit() else (1, s)


def _composition_batches(total: int, k: int, batch_size: int = _BATCH):
    """以 numpy 分块产出把 total 划分为 k 份正整数的全部组合。"""
    if k == 1:
        yield np.array([[total]], dtype=np.int64)
        return
    if k == 2:
        a = np.arange(1, total, dtype=np.int64)
        yield np.column_stack([a, total - a])
        return
    if k == 3:
        buf: List[np.ndarray] = []
        size = 0
        for a in range(1, total - 1):
            m = total - a - 1
            b = np.arange(1, m + 1, dtype=np.int64)
            part = np.empty((m, 3), dtype=np.int64)
            part[:, 0] = a
            part[:, 1] = b
            part[:, 2] = total - a - b
            buf.append(part)
            size += m
            if size >= batch_size:
                yield np.concatenate(buf)
                buf, size = [], 0
        if buf:
            yield np.concatenate(buf)
        return
    raise ValueError(f"不支持的供样人数: {k}")


def _features_for_subset(
    subset: Sequence[str],
    genotypes: Dict[str, Dict[str, Tuple[str, str]]],
    loci: Sequence[str],
    observed: Dict[str, Dict[str, int]],
) -> List[Tuple[str, str]]:
    """该子集需要核对的 (位点, 等位基因) 特征：观测等位基因 ∪ 子集成员基因型等位基因。"""
    feats: List[Tuple[str, str]] = []
    for locus in loci:
        alleles = set(observed[locus])
        for cid in subset:
            alleles.update(genotypes[cid][locus])
        for allele in sorted(alleles, key=_natural_key):
            feats.append((locus, allele))
    return feats


def _dosage_matrix(subset, genotypes, feats) -> np.ndarray:
    """G[i, j] = 子集第 i 人在特征 j 上的等位基因剂量（0/1/2）。"""
    G = np.zeros((len(subset), len(feats)), dtype=np.int64)
    for i, cid in enumerate(subset):
        for j, (locus, allele) in enumerate(feats):
            pair = genotypes[cid][locus]
            G[i, j] = (pair[0] == allele) + (pair[1] == allele)
    return G


def _describe(
    vector: Vector,
    genotypes: Dict[str, Dict[str, Tuple[str, str]]],
    loci: Sequence[str],
    observed: Dict[str, Dict[str, int]],
) -> dict:
    """把一份最优向量展开为可复核的逐位点逐等位基因明细。"""
    per_locus = []
    total2 = 0
    for locus in loci:
        alleles = set(observed[locus])
        for cid, _ in vector:
            alleles.update(genotypes[cid][locus])
        rows = []
        for allele in sorted(alleles, key=_natural_key):
            obs = observed[locus].get(allele, 0)
            pred2 = sum(
                copies
                * ((genotypes[cid][locus][0] == allele) + (genotypes[cid][locus][1] == allele))
                for cid, copies in vector
            )
            r2 = pred2 - 2 * obs
            total2 += abs(r2)
            rows.append(
                {
                    "allele": allele,
                    "observed": obs,
                    "predicted": pred2 / 2,
                    "residual": r2 / 2,
                }
            )
        per_locus.append({"locus": locus, "alleles": rows})
    return {
        "contributors": [{"candidateId": cid, "copies": c} for cid, c in vector],
        "residualSum": total2 / 2,
        "loci": per_locus,
    }


def solve(
    candidate_ids: Sequence[str],
    genotypes: Dict[str, Dict[str, Tuple[str, str]]],
    loci: Sequence[str],
    observed: Dict[str, Dict[str, int]],
    tolerance: int,
) -> dict:
    """在完整候选人子集 × 正整数份数的联合空间中裁决。

    返回 dict，status ∈ {"unique", "ambiguous", "no_solution"}。
    """
    depth = sum(observed[loci[0]].values())
    tol2 = 2 * tolerance

    best_k = None
    best_sum2 = 0
    best_count = 0
    best_witnesses: List[Vector] = []

    ordered_ids = sorted(candidate_ids, key=_natural_key)

    for k in (1, 2, 3):
        if depth < k:
            continue
        k_sum2 = None
        k_count = 0
        k_witnesses: List[Vector] = []

        for subset in combinations(ordered_ids, k):
            feats = _features_for_subset(subset, genotypes, loci, observed)
            o2 = np.array(
                [2 * observed[locus].get(allele, 0) for locus, allele in feats],
                dtype=np.int64,
            )
            G = _dosage_matrix(subset, genotypes, feats)

            for comps in _composition_batches(depth, k):
                pred = comps @ G                       # (n, 特征数)，双倍预测
                resid2 = np.abs(pred - o2)             # 双倍残差
                sums2 = resid2.sum(axis=1)
                feasible = resid2.max(axis=1) <= tol2  # 每个等位基因都须在上限内
                if not feasible.any():
                    continue
                fidx = np.nonzero(feasible)[0]
                fsums = sums2[fidx]
                local_min = int(fsums.min())
                if k_sum2 is not None and local_min > k_sum2:
                    continue
                local_count = int((fsums == local_min).sum())
                if k_sum2 is None or local_min < k_sum2:
                    k_sum2 = local_min
                    k_count = 0
                    k_witnesses = []
                k_count += local_count
                if len(k_witnesses) < 2:
                    take = np.nonzero(fsums == local_min)[0][: 2 - len(k_witnesses)]
                    for t in take:
                        comp = comps[fidx[t]]
                        vec: Vector = tuple(
                            sorted(
                                zip(subset, (int(x) for x in comp)),
                                key=lambda p: _natural_key(p[0]),
                            )
                        )
                        k_witnesses.append(vec)

        if k_sum2 is not None:
            # 供样人数最少者优先：一旦某 k 可行，不再尝试更大的 k
            best_k = k
            best_sum2 = k_sum2
            best_count = k_count
            best_witnesses = k_witnesses
            break

    base = {"totalDepth": depth, "tolerance": tolerance}
    if best_k is None:
        return {
            **base,
            "status": "no_solution",
            "objective": None,
            "optimalVectorCount": 0,
            "solution": None,
            "witnesses": [],
        }

    objective = {"contributors": best_k, "residualSum": best_sum2 / 2}
    if best_count == 1:
        return {
            **base,
            "status": "unique",
            "objective": objective,
            "optimalVectorCount": 1,
            "solution": _describe(best_witnesses[0], genotypes, loci, observed),
            "witnesses": None,
        }
    return {
        **base,
        "status": "ambiguous",
        "objective": objective,
        "optimalVectorCount": best_count,
        "solution": None,
        "witnesses": [_describe(v, genotypes, loci, observed) for v in best_witnesses[:2]],
    }
