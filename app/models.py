"""POST /api/mixture-adjudications 请求体校验。

非法输入（人数/位点数越界、编号重复、基因型缺失、各位点观测总深度不一致、
深度越界等）一律在此以 422 拒绝，绝不进入求解器、更不会伪装成唯一解。
"""

from typing import Dict, List, Tuple

from pydantic import BaseModel, Field, field_validator, model_validator

from .solver import MAX_DEPTH, MAX_OBSERVED_ALLELES

MAX_ALLELE_LEN = 24
MAX_NAME_LEN = 32
MAX_COUNT = 100_000


def _clean_name(value, what: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{what}必须是字符串")
    v = value.strip()
    if not v:
        raise ValueError(f"{what}不能为空")
    if len(v) > MAX_NAME_LEN:
        raise ValueError(f"{what}长度不能超过 {MAX_NAME_LEN} 个字符")
    return v


def _clean_allele(value) -> str:
    if not isinstance(value, str):
        raise ValueError("等位基因名称必须是字符串")
    v = value.strip()
    if not v:
        raise ValueError("等位基因名称不能为空")
    if len(v) > MAX_ALLELE_LEN:
        raise ValueError(f"等位基因名称长度不能超过 {MAX_ALLELE_LEN} 个字符")
    return v


class Candidate(BaseModel):
    id: str
    genotypes: Dict[str, Tuple[str, str]]  # 位点 -> 恰好两条等位基因

    @field_validator("id")
    @classmethod
    def _id_ok(cls, v):
        return _clean_name(v, "候选人编号")

    @field_validator("genotypes")
    @classmethod
    def _genotypes_ok(cls, g):
        cleaned = {}
        for locus, pair in g.items():
            name = _clean_name(locus, "基因型位点")
            if name in cleaned:
                raise ValueError(f"基因型位点重复: {name}")
            cleaned[name] = (_clean_allele(pair[0]), _clean_allele(pair[1]))
        return cleaned


class AdjudicationRequest(BaseModel):
    candidates: List[Candidate] = Field(min_length=4, max_length=8)
    loci: List[str] = Field(min_length=3, max_length=6)
    observed: Dict[str, Dict[str, int]]  # 位点 -> 等位基因 -> 观测拷贝数
    tolerance: int = Field(ge=0, le=MAX_DEPTH)

    @field_validator("loci")
    @classmethod
    def _loci_ok(cls, v):
        names = [_clean_name(x, "位点名称") for x in v]
        if len(set(names)) != len(names):
            raise ValueError("位点名称必须唯一")
        return names

    @field_validator("observed")
    @classmethod
    def _observed_ok(cls, obs):
        if not isinstance(obs, dict):
            raise ValueError("观测数据必须是按位点组织的对象")
        cleaned = {}
        for locus, table in obs.items():
            name = _clean_name(locus, "观测位点")
            if name in cleaned:
                raise ValueError(f"观测位点重复: {name}")
            if not isinstance(table, dict) or not table:
                raise ValueError(f"位点 {name} 至少需要一个等位基因观测")
            if len(table) > MAX_OBSERVED_ALLELES:
                raise ValueError(
                    f"位点 {name} 的等位基因种数超过上限 {MAX_OBSERVED_ALLELES}"
                )
            t = {}
            for allele, count in table.items():
                a = _clean_allele(allele)
                if a in t:
                    raise ValueError(f"位点 {name} 的等位基因重复: {a}")
                if isinstance(count, bool) or not isinstance(count, int):
                    raise ValueError(
                        f"位点 {name} 等位基因 {a} 的观测拷贝数必须是整数"
                    )
                if not 0 <= count <= MAX_COUNT:
                    raise ValueError(
                        f"位点 {name} 等位基因 {a} 的观测拷贝数须在 0~{MAX_COUNT} 之间"
                    )
                t[a] = count
            cleaned[name] = t
        return cleaned

    @model_validator(mode="after")
    def _cross_check(self):
        ids = [c.id for c in self.candidates]
        if len(set(ids)) != len(ids):
            raise ValueError("候选人编号必须唯一")
        locus_set = set(self.loci)
        for c in self.candidates:
            if set(c.genotypes) != locus_set:
                raise ValueError(f"候选人 {c.id} 的基因型位点与请求位点集合不一致")
        if set(self.observed) != locus_set:
            raise ValueError("观测数据必须覆盖且仅覆盖全部请求位点")
        totals = {locus: sum(self.observed[locus].values()) for locus in self.loci}
        if len(set(totals.values())) != 1:
            raise ValueError("所有位点的观测总深度必须一致")
        depth = next(iter(totals.values()))
        if not 1 <= depth <= MAX_DEPTH:
            raise ValueError(f"观测总深度必须在 1~{MAX_DEPTH} 之间")
        return self
