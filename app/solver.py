"""混合 DNA 样本裁决求解器。

模型约定（供复核）：
- 每位被选供样人分配正整数份数 n_i，份数之和等于各位点共同的观测总深度 D。
- 在某位点上，供样人的每份模板按 1/2 比例分摊到其两条等位基因
  （纯合子两条相同，即该等位基因获得 1 份）。因此等位基因 a 的预测拷贝数为
      pred(a) = Σ_i n_i * m_i(a) / 2,  m_i(a) ∈ {0, 1, 2}
  从而每位点预测总量 = Σ_i n_i = D，与观测总深度一致。
- 可行性：每个位点每种等位基因（观测与候选基因型的并集）都满足
      |pred(a) - obs(a)| <= tolerance
- 目标：先最少化供样人数 k ∈ {1, 2, 3}，再最小化全部位点全部等位基因的
  绝对残差和 Σ|pred - obs|。
- 唯一性：在最优目标值下，若按候选人编号排序的份数向量只有一种则为唯一；
  否则为歧义，须给出两份不同见证。

内部统一按“二倍拷贝数”整数计算（pred2 = Σ n_i*m_i, obs2 = 2*obs,
tol2 = 2*tolerance），避免半份带来的浮点误差，展示时再除以 2。
"""

from __future__ import annotations

from itertools import combinations

import numpy as np

# 单次裁决支持的观测总深度上限（保证组合枚举在有限时间内完成）。
MAX_TOTAL_DEPTH = 500
# 每个位点参与比较的等位基因数量上限（观测与候选基因型的并集）。
MAX_ALLELES_PER_LOCUS = 40

_BIG = np.int64(1 << 60)


def solve(candidates: list[dict], loci: list[dict], tolerance: int) -> dict:
    """在完整候选人子集与正整数份数的联合空间中执行裁决。

    candidates: [{"id": str, "genotypes": {locus_name: [a1, a2]}}]
    loci:       [{"name": str, "observed": {allele: count}}]
    返回可 JSON 序列化的裁决结果字典。
    """
    ordered = sorted(candidates, key=lambda c: c["id"])
    id_to_idx = {c["id"]: i for i, c in enumerate(ordered)}

    # 每个位点：等位基因并集（观测 ∪ 全部候选人基因型）与二倍观测向量。
    locus_blocks = []  # (locus_name, [allele, ...])
    obs2_parts = []
    for locus in loci:
        allele_set = set(locus["observed"])
        for c in ordered:
            allele_set.update(c["genotypes"][locus["name"]])
        alleles = sorted(allele_set)
        locus_blocks.append((locus["name"], alleles))
        obs2_parts.append(
            np.array([2 * locus["observed"].get(a, 0) for a in alleles], dtype=np.int64)
        )
    obs2_flat = np.concatenate(obs2_parts)
    depth = int(sum(loci[0]["observed"].values()))

    # 每位候选人的二倍等位基因计数向量（拼接全部位点，取值 0/1/2）。
    cand2 = []
    for c in ordered:
        parts = []
        for name, alleles in locus_blocks:
            index = {a: i for i, a in enumerate(alleles)}
            v = np.zeros(len(alleles), dtype=np.int64)
            a1, a2 = c["genotypes"][name]
            v[index[a1]] += 1
            v[index[a2]] += 1
            parts.append(v)
        cand2.append(np.concatenate(parts))

    tol2 = 2 * tolerance
    subsets_examined = 0
    assignments_examined = 0

    def record(level: list, vec: tuple, residual2: int) -> None:
        """level = [最优残差2 或 None, 见证向量列表(至多 2 个)]。"""
        if level[0] is None or residual2 < level[0]:
            level[0] = residual2
            level[1] = [vec]
        elif residual2 == level[0] and len(level[1]) < 2 and vec not in level[1]:
            level[1].append(vec)

    def needed(level: list, residual2: int) -> int:
        """当前残差下还需要收集的见证数量（0 表示可跳过）。"""
        if level[0] is None or residual2 < level[0]:
            return 2
        if residual2 == level[0]:
            return max(0, 2 - len(level[1]))
        return 0

    best = None  # (k, residual2)
    best_vectors: list[tuple] = []

    for k in (1, 2, 3):
        if depth < k:
            break  # 正整数份数要求份数和 >= k
        level: list = [None, []]
        for subset in combinations(range(len(ordered)), k):
            subsets_examined += 1
            mats = [cand2[i] for i in subset]
            ids = [ordered[i]["id"] for i in subset]

            if k == 1:
                assignments_examined += 1
                diff = np.abs(mats[0] * depth - obs2_flat)
                if int(diff.max()) <= tol2:
                    record(level, ((ids[0], depth),), int(diff.sum()))

            elif k == 2:
                x = np.arange(1, depth, dtype=np.int64)
                assignments_examined += int(x.size)
                y = depth - x
                pred = x[:, None] * mats[0] + y[:, None] * mats[1]
                diff = np.abs(pred - obs2_flat)
                feasible = diff.max(axis=1) <= tol2
                if not feasible.any():
                    continue
                masked = np.where(feasible, diff.sum(axis=1), _BIG)
                rmin = int(masked.min())
                need = needed(level, rmin)
                if need:
                    for pos in np.flatnonzero(masked == rmin)[:need]:
                        record(
                            level,
                            ((ids[0], int(x[pos])), (ids[1], int(y[pos]))),
                            rmin,
                        )

            else:  # k == 3
                m0, m1, m2 = mats
                for xv in range(1, depth - 1):
                    ys = np.arange(1, depth - xv, dtype=np.int64)
                    assignments_examined += int(ys.size)
                    zs = depth - xv - ys
                    pred = xv * m0 + ys[:, None] * m1 + zs[:, None] * m2
                    diff = np.abs(pred - obs2_flat)
                    feasible = diff.max(axis=1) <= tol2
                    if not feasible.any():
                        continue
                    masked = np.where(feasible, diff.sum(axis=1), _BIG)
                    rmin = int(masked.min())
                    need = needed(level, rmin)
                    if not need:
                        continue
                    for pos in np.flatnonzero(masked == rmin)[:need]:
                        record(
                            level,
                            (
                                (ids[0], xv),
                                (ids[1], int(ys[pos])),
                                (ids[2], int(zs[pos])),
                            ),
                            rmin,
                        )

        if level[0] is not None:
            best = (k, level[0])
            best_vectors = level[1]
            break  # 供样人数已最少，不再搜索更大子集

    result = {
        "status": "no_solution",
        "totalDepth": depth,
        "tolerance": tolerance,
        "objective": None,
        "solution": None,
        "witnesses": None,
        "examined": {
            "subsets": subsets_examined,
            "assignments": assignments_examined,
        },
    }
    if best is None:
        return result

    k, residual2 = best
    result["objective"] = {
        "contributors": k,
        "absoluteResidual": residual2 / 2,
    }
    witnesses = [
        _witness(vec, id_to_idx, cand2, locus_blocks, obs2_parts)
        for vec in best_vectors
    ]
    if len(witnesses) == 1:
        result["status"] = "unique"
        result["solution"] = witnesses[0]
    else:
        result["status"] = "ambiguous"
        result["witnesses"] = witnesses
    return result


def _witness(
    vector: tuple,
    id_to_idx: dict,
    cand2: list,
    locus_blocks: list,
    obs2_parts: list,
) -> dict:
    """把一份最优份数向量展开为可复核的逐位点预测/观测对照。"""
    contributors = [{"id": cid, "copies": n} for cid, n in vector]
    per_locus = []
    total2 = 0
    pos = 0
    for (name, alleles), obs2 in zip(locus_blocks, obs2_parts):
        width = len(alleles)
        pred2 = np.zeros(width, dtype=np.int64)
        for cid, n in vector:
            pred2 += n * cand2[id_to_idx[cid]][pos : pos + width]
        rows = []
        residual2 = 0
        for allele, p2, o2 in zip(alleles, pred2, obs2):
            d2 = abs(int(p2) - int(o2))
            residual2 += d2
            rows.append(
                {
                    "allele": allele,
                    "observed": int(o2) // 2,
                    "predicted": int(p2) / 2,
                    "absDiff": d2 / 2,
                }
            )
        total2 += residual2
        per_locus.append(
            {"locus": name, "residual": residual2 / 2, "alleles": rows}
        )
        pos += width
    return {
        "contributors": contributors,
        "absoluteResidual": total2 / 2,
        "perLocus": per_locus,
    }
