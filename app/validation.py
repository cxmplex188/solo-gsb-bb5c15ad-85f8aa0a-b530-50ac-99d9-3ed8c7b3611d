"""跨字段语义校验：位点越界、重复位点、酶名称重复等。

校验失败统一抛出 :class:`DigestValidationError`，由主应用转换为
统一结构的 HTTP 422 响应（与 Pydantic 字段校验同一响应格式）。
"""
from fastapi import HTTPException

from app.config import settings
from app.models import DigestRequest


class DigestValidationError(HTTPException):
    def __init__(self, errors: list[dict]):
        super().__init__(status_code=422, detail={"code": "validation_failed", "errors": errors})


def _err(loc: list[str], message: str, code: str = "invalid") -> dict:
    return {"loc": loc, "msg": message, "type": code}


def validate_request(req: DigestRequest) -> DigestRequest:
    errors: list[dict] = []
    length = req.plasmid_length

    if req.min_visible_fragment >= length:
        errors.append(_err(
            ["min_visible_fragment"],
            f"最小可见片段长度（{req.min_visible_fragment}）必须小于质粒长度（{length}）",
            "min_visible_too_large",
        ))

    if len(req.enzymes) > settings.MAX_ENZYMES:
        errors.append(_err(
            ["enzymes"],
            f"酶数量（{len(req.enzymes)}）超过允许上限 {settings.MAX_ENZYMES}",
            "too_many_enzymes",
        ))

    seen_names: set[str] = set()
    for i, enzyme in enumerate(req.enzymes):
        name = enzyme.name
        if name in seen_names:
            errors.append(_err(
                ["enzymes", i, "name"],
                f"酶名称重复：{name}",
                "duplicate_enzyme_name",
            ))
        seen_names.add(name)

        for field in ("reference_sites", "variant_sites"):
            sites = getattr(enzyme, field)
            sample = "reference" if field == "reference_sites" else "variant"

            out_of_range = [s for s in sites if s < 0 or s >= length]
            if out_of_range:
                errors.append(_err(
                    ["enzymes", i, field],
                    f"酶 {name} 在{('参考型' if sample == 'reference' else '变异型')}"
                    f"中存在越界切点 {out_of_range}：合法范围为 0 ~ {length - 1}",
                    "site_out_of_range",
                ))

            unique_sites = set(sites)
            if len(unique_sites) != len(sites):
                dupes = sorted({s for s in sites if sites.count(s) > 1})
                errors.append(_err(
                    ["enzymes", i, field],
                    f"酶 {name} 在{('参考型' if sample == 'reference' else '变异型')}"
                    f"中存在重复切点 {dupes}",
                    "duplicate_site",
                ))

    if errors:
        raise DigestValidationError(errors)
    return req
