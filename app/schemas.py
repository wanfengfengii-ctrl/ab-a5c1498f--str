"""请求体的结构校验。非法输入一律以 422 拒绝，绝不进入求解器。"""

from __future__ import annotations

from pydantic import BaseModel, Field, field_validator, model_validator

from .solver import MAX_ALLELES_PER_LOCUS, MAX_TOTAL_DEPTH

LABEL_MAX = 24
MIN_CANDIDATES = 4
MAX_CANDIDATES = 8
MIN_LOCI = 3
MAX_LOCI = 6


def _clean_label(value: object, what: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{what}必须是字符串")
    text = value.strip()
    if not text:
        raise ValueError(f"{what}不能为空")
    if len(text) > LABEL_MAX:
        raise ValueError(f"{what}长度不能超过 {LABEL_MAX} 个字符")
    return text


class LocusInput(BaseModel):
    name: str
    observed: dict[str, int]

    @field_validator("name")
    @classmethod
    def _name_ok(cls, v: str) -> str:
        return _clean_label(v, "位点名称")

    @field_validator("observed")
    @classmethod
    def _observed_ok(cls, v: dict[str, int]) -> dict[str, int]:
        if not v:
            raise ValueError("每个位点至少需要一种观测等位基因")
        if len(v) > MAX_ALLELES_PER_LOCUS:
            raise ValueError(f"每个位点的观测等位基因数不能超过 {MAX_ALLELES_PER_LOCUS}")
        cleaned: dict[str, int] = {}
        for allele, count in v.items():
            name = _clean_label(allele, "等位基因")
            if isinstance(count, bool) or not isinstance(count, int):
                raise ValueError(f"等位基因 {name} 的观测拷贝数必须为整数")
            if count < 0:
                raise ValueError(f"等位基因 {name} 的观测拷贝数不能为负数")
            if count > MAX_TOTAL_DEPTH:
                raise ValueError(
                    f"等位基因 {name} 的观测拷贝数超过上限 {MAX_TOTAL_DEPTH}"
                )
            cleaned[name] = count
        return cleaned


class CandidateInput(BaseModel):
    id: str
    genotypes: dict[str, list[str]]

    @field_validator("id")
    @classmethod
    def _id_ok(cls, v: str) -> str:
        return _clean_label(v, "候选人编号")

    @field_validator("genotypes")
    @classmethod
    def _genotypes_ok(cls, v: dict[str, list[str]]) -> dict[str, list[str]]:
        if not v:
            raise ValueError("候选人必须提供基因型")
        cleaned: dict[str, list[str]] = {}
        for locus, alleles in v.items():
            name = _clean_label(locus, "位点名称")
            if not isinstance(alleles, list) or len(alleles) != 2:
                raise ValueError(f"位点 {name} 必须提供恰好两条等位基因")
            cleaned[name] = [_clean_label(a, "等位基因") for a in alleles]
        return cleaned


class AdjudicationRequest(BaseModel):
    tolerance: int = Field(ge=0, le=1_000_000)
    candidates: list[CandidateInput] = Field(
        min_length=MIN_CANDIDATES, max_length=MAX_CANDIDATES
    )
    loci: list[LocusInput] = Field(min_length=MIN_LOCI, max_length=MAX_LOCI)

    @model_validator(mode="after")
    def _cross_check(self) -> "AdjudicationRequest":
        ids = [c.id for c in self.candidates]
        if len(set(ids)) != len(ids):
            raise ValueError("候选人编号必须唯一")
        names = [l.name for l in self.loci]
        if len(set(names)) != len(names):
            raise ValueError("位点名称必须唯一")
        locus_set = set(names)
        for c in self.candidates:
            if set(c.genotypes) != locus_set:
                raise ValueError(
                    f"候选人 {c.id} 的基因型位点集合必须与位点列表完全一致"
                )
        depths = {sum(l.observed.values()) for l in self.loci}
        if len(depths) != 1:
            raise ValueError("各位点的观测总深度不一致，无法裁决")
        depth = depths.pop()
        if depth < 1:
            raise ValueError("各位点的观测总深度必须为正整数")
        if depth > MAX_TOTAL_DEPTH:
            raise ValueError(f"观测总深度 {depth} 超过上限 {MAX_TOTAL_DEPTH}")
        for locus in self.loci:
            allele_set = set(locus.observed)
            for c in self.candidates:
                allele_set.update(c.genotypes[locus.name])
            if len(allele_set) > MAX_ALLELES_PER_LOCUS:
                raise ValueError(
                    f"位点 {locus.name} 的等位基因并集超过上限 {MAX_ALLELES_PER_LOCUS}"
                )
        return self
