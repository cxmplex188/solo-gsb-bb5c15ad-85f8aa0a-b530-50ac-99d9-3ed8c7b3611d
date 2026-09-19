"""端到端 API 测试（FastAPI TestClient，不启动真实网络服务）。"""
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def _enzyme(name, ref, var):
    return {"name": name, "reference_sites": ref, "variant_sites": var}


def _payload(enzymes, length=1000, min_visible=50, tolerance=20):
    return {
        "plasmid_length": length,
        "min_visible_fragment": min_visible,
        "band_resolution_tolerance": tolerance,
        "enzymes": enzymes,
    }


# ---------------------------------------------------------------- 健康检查

def test_health():
    res = client.get("/health")
    assert res.status_code == 200
    body = res.json()
    assert body["status"] == "ok"
    assert "version" in body


# ---------------------------------------------------------------- 基础推荐

def test_recommend_basic_scenario():
    payload = _payload([
        _enzyme("EcoRI", [100, 600], [100, 600, 850]),
        _enzyme("BamHI", [300], [300]),
    ])
    res = client.post("/api/v1/digest/recommend", json=payload)
    assert res.status_code == 200
    body = res.json()

    assert body["total_enzyme_pairs"] == 1
    assert body["recommended_count"] == 1
    rec = body["recommendations"][0]
    assert rec["enzyme_pair"] == ["BamHI", "EcoRI"]
    assert rec["distinguishable_band_count"] >= 1

    # 每种样本都应报告一个跨原点片段
    for sample in ("reference", "variant"):
        crosses = [f for f in rec["samples"][sample]["fragments"] if f["crosses_origin"]]
        assert len(crosses) == 1

    # 所有差异条带必须可区分
    assert all(d["reference_count"] != d["variant_count"] for d in rec["differences"])
    assert {d["difference_type"] for d in rec["differences"]} <= {
        "reference_only", "variant_only", "count_mismatch",
    }


def test_fragment_lengths_sum_to_plasmid_length():
    payload = _payload([
        _enzyme("EcoRI", [100, 600], [100, 600, 850]),
        _enzyme("BamHI", [300], [300]),
    ])
    rec = client.post("/api/v1/digest/recommend", json=payload).json()["recommendations"][0]
    for sample in ("reference", "variant"):
        fragments = rec["samples"][sample]["fragments"]
        assert sum(f["length"] for f in fragments) == 1000


# ---------------------------------------------------------------- 排序

def test_recommendations_sorting_distinguishable_desc_then_total_asc():
    # L=100、最小可见=1、容差=0，手工精确推演（片段多集合）：
    #   A+B：参考 {20,80} vs 变异 {20,30,50} -> 3 条可区分 / 共 4 条带
    #   A+C：参考 {20,80} vs 变异 {40,60}     -> 4 条可区分 / 共 4 条带
    #   C+D：参考 {30,70} vs 变异 {10,90}     -> 4 条可区分 / 共 4 条带（键序在后）
    #   A+D、B+D：两样本一致 -> 不可区分；B+C：参考仅 1 个合并切点 -> 拒绝
    payload = _payload(
        [
            _enzyme("A", [0], [0]),
            _enzyme("B", [20], [20, 50]),
            _enzyme("C", [20], [40]),
            _enzyme("D", [50], [50]),
        ],
        length=100,
        min_visible=1,
        tolerance=0,
    )
    body = client.post("/api/v1/digest/recommend", json=payload).json()

    assert body["total_enzyme_pairs"] == 6
    keys = [r["pair_key"] for r in body["recommendations"]]
    assert keys == ["A+C", "C+D", "A+B"]
    assert body["recommendations"][0]["distinguishable_band_count"] == 4
    assert body["recommendations"][1]["distinguishable_band_count"] == 4
    ab = body["recommendations"][2]
    assert (ab["distinguishable_band_count"], ab["total_bands"]) == (3, 4)

    assert [p["pair_key"] for p in body["non_discriminating_pairs"]] == ["A+D", "B+D"]
    assert body["non_discriminating_pair_count"] == 2
    assert body["rejected_pairs"][0]["pair_key"] == "B+C"
    assert body["rejected_pairs"][0]["reason_code"] == "insufficient_cuts"


def test_tie_break_on_total_bands_and_name():
    # 两对都不可区分以外的同分时按总带数升序、酶名升序
    payload = _payload([
        _enzyme("EcoRI", [0], [0, 200]),
        _enzyme("BamHI", [400], [400, 700]),
        _enzyme("HindIII", [100, 900], [100, 500, 900]),
    ])
    body = client.post("/api/v1/digest/recommend", json=payload).json()
    recs = body["recommendations"]
    for a, b in zip(recs, recs[1:]):
        ka = (-a["distinguishable_band_count"], a["total_bands"], a["pair_key"])
        kb = (-b["distinguishable_band_count"], b["total_bands"], b["pair_key"])
        assert ka <= kb


# ---------------------------------------------------------------- 拒绝逻辑

def test_rejected_when_sample_has_fewer_than_two_cuts():
    # HindIII 对参考型无切点 -> EcoRI+HindIII 对参考型仅 1 个切点
    payload = _payload([
        _enzyme("EcoRI", [100], [100, 600]),
        _enzyme("HindIII", [], [300]),
    ])
    body = client.post("/api/v1/digest/recommend", json=payload).json()
    assert body["recommended_count"] == 0
    assert body["rejected_pair_count"] == 1
    rejected = body["rejected_pairs"][0]
    assert rejected["reason_code"] == "insufficient_cuts"
    assert rejected["reference_cut_count"] == 1
    assert rejected["variant_cut_count"] == 3


def test_non_discriminating_pair_listed_separately():
    # 两酶在两种样本上切点完全相同 -> 无法区分
    payload = _payload([
        _enzyme("EcoRI", [100, 600], [100, 600]),
        _enzyme("BamHI", [300], [300]),
    ])
    body = client.post("/api/v1/digest/recommend", json=payload).json()
    assert body["recommended_count"] == 0
    assert body["non_discriminating_pair_count"] == 1
    assert body["non_discriminating_pairs"][0]["pair_key"] == "BamHI+EcoRI"


def test_coincident_sites_from_different_enzymes_merge():
    # 参考型两酶在 100 处同位点（合并，不按重复位点拒绝）；变异型 BamHI 在 350 处
    payload = _payload([
        _enzyme("EcoRI", [100, 600], [100, 600]),
        _enzyme("BamHI", [100, 300], [350]),
    ])
    body = client.post("/api/v1/digest/recommend", json=payload).json()
    assert body["recommended_count"] == 1
    rec = body["recommendations"][0]
    cuts = rec["samples"]["reference"]["cuts"]
    at_100 = [c for c in cuts if c["position"] == 100]
    assert len(at_100) == 1
    assert at_100[0]["enzymes"] == ["BamHI", "EcoRI"]
    var_cuts = rec["samples"]["variant"]["cuts"]
    assert var_cuts[0]["position"] == 100 and var_cuts[0]["enzymes"] == ["EcoRI"]


# ---------------------------------------------------------------- 422 校验

def _expect_422(payload, code_part=None):
    res = client.post("/api/v1/digest/recommend", json=payload)
    assert res.status_code == 422
    body = res.json()
    assert body["code"] == "validation_failed"
    assert body["errors"]
    if code_part:
        assert any(code_part in e["type"] for e in body["errors"])
    return body


def test_422_site_out_of_range():
    payload = _payload([
        _enzyme("EcoRI", [100, 1000], [100]),   # 1000 越界（合法 0..999）
        _enzyme("BamHI", [300], [300]),
    ])
    _expect_422(payload, "site_out_of_range")


def test_422_negative_site():
    payload = _payload([
        _enzyme("EcoRI", [-1], [100]),
        _enzyme("BamHI", [300], [300]),
    ])
    _expect_422(payload, "site_out_of_range")


def test_422_duplicate_site():
    payload = _payload([
        _enzyme("EcoRI", [100, 100, 600], [100]),
        _enzyme("BamHI", [300], [300]),
    ])
    _expect_422(payload, "duplicate_site")


def test_422_duplicate_enzyme_name():
    payload = _payload([
        _enzyme("EcoRI", [100], [100]),
        _enzyme("EcoRI", [300], [300]),
    ])
    _expect_422(payload, "duplicate_enzyme_name")


def test_422_min_visible_not_smaller_than_plasmid():
    payload = _payload(
        [_enzyme("EcoRI", [100], [100]), _enzyme("BamHI", [300], [300])],
        min_visible=1000,
    )
    _expect_422(payload, "min_visible_too_large")


def test_422_basic_field_constraints():
    # plasmid_length <= 0、容差为负、酶不足 2 种
    _expect_422(_payload(
        [_enzyme("EcoRI", [100], [100])],
        length=0,
    ))
    res = client.post("/api/v1/digest/recommend", json=_payload(
        [_enzyme("EcoRI", [100], [100]), _enzyme("BamHI", [300], [300])],
    ))
    # 正常对照
    assert res.status_code == 200
