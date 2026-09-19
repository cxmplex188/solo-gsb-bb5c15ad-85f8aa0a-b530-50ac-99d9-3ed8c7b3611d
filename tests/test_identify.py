"""单元测试：预测条带与观测条带的容差内一对一最小代价匹配。"""
import itertools
import math

from app.identify import identify, match_sample
from app.models import IdentifyRequest


# ---------------------------------------------------------------- 基本匹配

def test_exact_and_near_matches_within_tolerance():
    result = match_sample([100, 300, 500], [102, 298, 505], tolerance=30, sample="reference")
    assert result.matched_count == 3
    assert result.missing_bands == []
    assert result.extra_bands == []
    assert result.total_error == 9.0
    assert [(m.predicted, m.observed) for m in result.matches] == [
        (100, 102), (300, 298), (500, 505),
    ]
    assert all(m.within_tolerance for m in result.matches)


def test_missing_and_extra_bands():
    result = match_sample([100, 300, 500], [100, 900], tolerance=20, sample="reference")
    assert result.matched_count == 1
    assert result.missing_band_count == 2
    assert result.extra_band_count == 1
    assert result.missing_bands == [300, 500]
    assert result.extra_bands == [900]
    assert result.total_error == 0.0


def test_empty_predicted_all_observed_extra():
    result = match_sample([], [100, 200], tolerance=20, sample="variant")
    assert result.matched_count == 0
    assert result.missing_bands == []
    assert result.extra_bands == [100, 200]
    assert result.total_error == 0.0


def test_empty_observed_all_predicted_missing():
    result = match_sample([100, 200], [], tolerance=20, sample="reference")
    assert result.matched_count == 0
    assert result.missing_bands == [100, 200]
    assert result.extra_bands == []


# ---------------------------------------------------------------- 最小代价

def test_matching_minimizes_total_cost_not_greedy():
    # 贪心按排序顺序会把 100 配 111（误差 11），110 失配；
    # 最优一对一为 100↔100、110↔111，总误差 1
    result = match_sample([100, 110], [100, 111], tolerance=20, sample="reference")
    assert [(m.predicted, m.observed) for m in result.matches] == [(100, 100), (110, 111)]
    assert result.total_error == 1.0


def test_maximum_cardinality_preferred_even_when_cost_higher():
    # 0 误差只能配 1 对（100↔100）；可配 2 对的方案误差 0+18=18，应优先配对数
    result = match_sample([100, 108], [100, 126], tolerance=20, sample="reference")
    assert result.matched_count == 2
    assert result.total_error == 18.0
    assert result.missing_band_count == 0
    assert result.extra_band_count == 0


def test_edge_at_exact_tolerance_boundary_is_matchable():
    result = match_sample([100], [110], tolerance=10, sample="reference")
    assert result.matched_count == 1
    assert result.total_error == 10.0

    strict = match_sample([100], [110.0001], tolerance=10, sample="reference")
    assert strict.matched_count == 0
    assert strict.missing_bands == [100]
    assert strict.extra_bands == [110.0001]


def test_input_order_does_not_affect_result():
    a = match_sample([500, 100, 300], [505, 102, 298], tolerance=30, sample="reference")
    b = match_sample([100, 300, 500], [298, 505, 102], tolerance=30, sample="reference")
    assert a.model_dump() == b.model_dump()


# ---------------------------------------------------------------- 暴力对照

def _brute_force_optimal(predicted, observed, tolerance):
    """枚举所有部分一对一配对：先最大配对数，再最小总绝对误差。"""
    n, m = len(predicted), len(observed)
    for r in range(min(n, m), -1, -1):
        best_cost = None
        best_pairs = None
        for ip in itertools.combinations(range(n), r):
            for jo in itertools.permutations(range(m), r):
                pairs = list(zip(ip, jo))
                if all(abs(predicted[i] - observed[j]) <= tolerance for i, j in pairs):
                    cost = sum(abs(predicted[i] - observed[j]) for i, j in pairs)
                    if best_cost is None or cost < best_cost:
                        best_cost = cost
                        best_pairs = pairs
        if best_pairs is not None:
            return r, best_cost, best_pairs
    return 0, 0.0, []


def test_matches_brute_force_on_random_cases():
    import random

    rng = random.Random(20260919)
    for _ in range(60):
        predicted = [rng.randint(0, 500) for _ in range(rng.randint(0, 4))]
        observed = [rng.randint(0, 500) for _ in range(rng.randint(0, 4))]
        tolerance = rng.choice([0, 10, 25, 60])
        result = match_sample(predicted, observed, tolerance, "reference")
        count, cost, _pairs = _brute_force_optimal(predicted, observed, tolerance)
        assert result.matched_count == count
        assert math.isclose(result.total_error, cost, abs_tol=1e-9)
        assert result.missing_band_count == len(predicted) - result.matched_count
        assert result.extra_band_count == len(observed) - result.matched_count


# ---------------------------------------------------------------- 判定

def test_identify_decision_reference_variant_and_ambiguous():
    req = IdentifyRequest(
        reference_bands=[100, 300, 500],
        variant_bands=[120, 320, 600],
        observed_bands=[102, 298, 505],
        tolerance=30,
    )
    resp = identify(req)
    assert resp.decision == "reference"
    assert resp.reference.total_error < resp.variant.total_error
    assert "参考型" in resp.reasoning

    req2 = req.model_copy(update={"observed_bands": [122, 318, 595]})
    assert identify(req2).decision == "variant"

    req3 = IdentifyRequest(
        reference_bands=[100], variant_bands=[102],
        observed_bands=[101], tolerance=5,
    )
    resp3 = identify(req3)
    assert resp3.decision == "ambiguous"
    assert resp3.reference.total_error == resp3.variant.total_error == 1.0
    assert "无法判定" in resp3.reasoning


def test_identify_uses_same_observed_list_for_both_samples():
    req = IdentifyRequest(
        reference_bands=[100, 300],
        variant_bands=[100, 300, 500],
        observed_bands=[100, 300, 500],
        tolerance=5,
    )
    resp = identify(req)
    # 变异型 3 条全部精确配对；参考型有 1 条额外条带 -> variant
    assert resp.variant.matched_count == 3
    assert resp.variant.extra_band_count == 0
    assert resp.reference.matched_count == 2
    assert resp.reference.extra_band_count == 1
    assert resp.decision == "variant"
