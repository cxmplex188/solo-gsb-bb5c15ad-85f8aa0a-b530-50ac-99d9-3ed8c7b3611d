"""请求 / 响应数据模型。

输入使用 Pydantic 做基础类型约束；跨字段的语义校验
（位点越界、重复位点、切点数量等）在 :mod:`app.validation` 中完成，
两者最终都输出统一结构的 422 响应。
"""
from typing import Annotated, Literal

from pydantic import BaseModel, Field

SampleName = Literal["reference", "variant"]

# 条带长度 / 匹配容差：非负有限浮点数（JSON 整数自动转换）
NonNegativeBand = Annotated[float, Field(ge=0, allow_inf_nan=False)]


# ---------------------------------------------------------------- 输入模型

class EnzymeInput(BaseModel):
    name: str = Field(..., min_length=1, description="酶名称，如 EcoRI")
    reference_sites: list[int] = Field(
        default_factory=list, description="参考型上该酶的候选切点位置（bp，0-based）"
    )
    variant_sites: list[int] = Field(
        default_factory=list, description="变异型上该酶的候选切点位置（bp，0-based）"
    )


class DigestRequest(BaseModel):
    plasmid_length: int = Field(..., gt=0, description="环状质粒总长度（bp）")
    min_visible_fragment: int = Field(
        ..., ge=0, description="最小可见片段长度（bp），短于此值的片段不参与条带聚类"
    )
    band_resolution_tolerance: int = Field(
        ..., ge=0, description="条带分辨容差（bp），长度差不超过该值的相邻片段聚为同一条带"
    )
    enzymes: list[EnzymeInput] = Field(..., min_length=2, description="参与枚举的酶，至少 2 种")


# ---------------------------------------------------------------- 输出模型

class CutPointOut(BaseModel):
    position: int
    enzymes: list[str] = Field(description="在该位置产生切点的酶（可能两酶同位点）")


class FragmentOut(BaseModel):
    start: int = Field(description="片段起始切点位置")
    end: int = Field(description="片段终止切点位置（跨原点时为环绕后的位置）")
    length: int
    crosses_origin: bool = Field(description="是否为跨过坐标原点（0 点）的片段")


class SampleDigestOut(BaseModel):
    cuts: list[CutPointOut]
    fragments: list[FragmentOut]
    visible_fragment_count: int
    hidden_fragment_count: int = Field(description="短于最小可见片段、未参与成带的片段数")


class BandMemberOut(FragmentOut):
    sample: SampleName


class BandOut(BaseModel):
    band_id: int = Field(description="条带编号，按代表长度升序，从 1 开始")
    representative_length: float = Field(description="条带内片段长度平均值（bp）")
    min_length: int
    max_length: int
    reference_count: int
    variant_count: int
    distinguishable: bool = Field(description="两种样本在该条带中的片段数是否不同")
    members: list[BandMemberOut]


class DifferenceOut(BaseModel):
    band_id: int
    representative_length: float
    difference_type: Literal["reference_only", "variant_only", "count_mismatch"]
    reference_count: int
    variant_count: int
    description: str = Field(description="中文差异依据说明")


class PairRecommendationOut(BaseModel):
    enzyme_pair: list[str]
    pair_key: str
    samples: dict[SampleName, SampleDigestOut]
    bands: list[BandOut]
    total_bands: int
    distinguishable_band_count: int
    differences: list[DifferenceOut]


class RejectedPairOut(BaseModel):
    enzyme_pair: list[str]
    pair_key: str
    reason_code: Literal["insufficient_cuts"]
    reason: str
    reference_cut_count: int
    variant_cut_count: int


class NonDiscriminatingPairOut(BaseModel):
    enzyme_pair: list[str]
    pair_key: str
    total_bands: int = Field(description="两种样本条带完全一致时的条带总数")


class DigestResponse(BaseModel):
    plasmid_length: int
    min_visible_fragment: int
    band_resolution_tolerance: int
    total_enzyme_pairs: int
    recommended_count: int
    rejected_pair_count: int
    non_discriminating_pair_count: int
    recommendations: list[PairRecommendationOut]
    rejected_pairs: list[RejectedPairOut]
    non_discriminating_pairs: list[NonDiscriminatingPairOut]


class HealthResponse(BaseModel):
    status: Literal["ok"]
    service: str
    version: str


# ---------------------------------------------------------------- 鉴别（identify）模型

class IdentifyRequest(BaseModel):
    reference_bands: list[NonNegativeBand] = Field(
        default_factory=list, description="参考型预测条带长度列表（bp）"
    )
    variant_bands: list[NonNegativeBand] = Field(
        default_factory=list, description="变异型预测条带长度列表（bp）"
    )
    observed_bands: list[NonNegativeBand] = Field(
        ..., min_length=1, description="一次凝胶观测到的条带长度列表（bp），至少 1 条"
    )
    band_match_tolerance: NonNegativeBand = Field(
        ..., description="条带匹配容差（bp）：长度差不超过该值的预测/观测条带才允许配对"
    )


class BandMatchOut(BaseModel):
    predicted_length: float = Field(description="预测条带长度（bp）")
    observed_length: float = Field(description="与之配对的观测条带长度（bp）")
    error: float = Field(description="预测与观测长度差的绝对值（bp）")
    description: str = Field(description="中文匹配依据说明")


class SampleMatchOut(BaseModel):
    predicted_band_count: int
    matched_count: int
    missing_count: int = Field(description="缺失条带数：未在观测中找到配对的预测条带")
    extra_count: int = Field(description="额外条带数：未匹配到任何预测条带的观测条带")
    total_error: float = Field(description="全部匹配对的长度差之和（bp）")
    matches: list[BandMatchOut] = Field(description="逐条匹配依据，按预测长度升序")
    missing_bands: list[float] = Field(description="缺失条带：未匹配到观测的预测长度，升序")
    extra_bands: list[float] = Field(description="额外条带：未匹配到预测的观测长度，升序")


class IdentifyResponse(BaseModel):
    band_match_tolerance: float
    observed_band_count: int
    samples: dict[SampleName, SampleMatchOut]
    verdict: Literal["reference", "variant", "ambiguous"]
    verdict_reason: str = Field(description="中文判定依据说明")
