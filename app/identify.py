"""凝胶条带鉴别：预测条带与一次凝胶观测条带的确定性一对一最小代价匹配。

匹配规则
--------
对每种样本（参考型 / 变异型）的预测条带长度列表 P 与**同一次**凝胶
观测列表 O，在给定容差 t 内做一对一匹配：

  * 仅当 ``|p - o| <= t`` 时，预测条带 p 与观测条带 o 之间才存在可行配对；
  * 每条预测、每条观测至多参与一次配对；
  * 先最大化可行配对数（匈牙利算法求最小权二分图匹配，未配对走代价为
    DUMMY 的虚拟节点，DUMMY 大于任一可行边权，故优先配对），再在配对数
    相同的方案中最小化总绝对误差；
  * 未被配对的预测条带记为"缺失条带"（凝胶上应有而未见），未被配对的
    观测条带记为"额外条带"（凝胶上出现但预测不到）。

两种样本各自完成匹配后，按与匹配目标一致的优先级比较：
**先比缺失+额外条带总数（容差内解释的条带数越多越优），总数相同时再比
总绝对误差；两者都相同为 ``ambiguous``**。仅比较配对误差会让"一条都没
配上"的样本以 0 误差与完全吻合的样本打平，产生退化判定。

确定性：实现为 O(k^3) 匈牙利算法，所有输入列表先按长度排序，匹配结果
只取决于输入数值，与列表原始顺序无关。
"""
import math

from app.models import (
    BandMatchOut,
    IdentifyRequest,
    IdentifyResponse,
    SampleIdentificationOut,
)

_SAMPLE_LABEL = {"reference": "参考型", "variant": "变异型"}


def _hungarian(cost: list[list[float]]) -> list[int]:
    """方阵最小权指派（1 基下标的 cp-algorithms 实现）。

    返回 ``assignment[i] = j``：第 i 行（1-based）指派给第 j 列（1-based）。
    """
    n = len(cost) - 1
    m = len(cost[0]) - 1
    u = [0.0] * (n + 1)
    v = [0.0] * (m + 1)
    p = [0] * (m + 1)   # p[j]：当前指派给第 j 列的行
    way = [0] * (m + 1)

    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = [math.inf] * (m + 1)
        used = [False] * (m + 1)
        while p[j0] != 0:
            used[j0] = True
            i0 = p[j0]
            delta = math.inf
            j1 = 0
            for j in range(1, m + 1):
                if not used[j]:
                    cur = cost[i0][j] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j] = cur
                        way[j] = j0
                    if minv[j] < delta:
                        delta = minv[j]
                        j1 = j
            for j in range(m + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break

    assignment = [0] * (n + 1)
    for j in range(1, m + 1):
        assignment[p[j]] = j
    return assignment


def match_sample(predicted: list[float], observed: list[float],
                 tolerance: float, sample: str) -> SampleIdentificationOut:
    """单种样本的预测条带与观测条带做容差内一对一最小代价匹配。"""
    preds = sorted(predicted)
    obs = sorted(observed)
    n, m = len(preds), len(obs)

    # 可行边（容差内）的权值之和上界；DUMMY 取 S+1，保证"多配一条可行边"
    # 始终比让该条带走未配对虚拟节点更优。
    feasible = [
        (i, j, abs(preds[i] - obs[j]))
        for i in range(n)
        for j in range(m)
        if abs(preds[i] - obs[j]) <= tolerance
    ]
    edge_sum = sum(w for _, _, w in feasible)
    dummy = edge_sum + 1.0

    # 扩展为 (n+m)×(n+m) 方阵：行=预测+虚拟行，列=观测+虚拟列
    k = n + m
    cost = [[0.0] * (k + 1) for _ in range(k + 1)]
    for i in range(1, k + 1):
        for j in range(1, k + 1):
            cost[i][j] = math.inf
    # 虚拟行/虚拟列之间自由配对
    for i in range(n + 1, k + 1):
        for j in range(m + 1, k + 1):
            cost[i][j] = 0.0
    # 真实预测 i 缺配 -> 它自己的虚拟列；真实观测 j 多配 -> 它自己的虚拟行
    for i in range(n):
        cost[i + 1][m + 1 + i] = dummy
    for j in range(m):
        cost[n + 1 + j][j + 1] = dummy
    # 容差内真实配对
    for i, j, w in feasible:
        cost[i + 1][j + 1] = w

    assignment = _hungarian(cost)

    matches: list[BandMatchOut] = []
    matched_pred: set[int] = set()
    matched_obs: set[int] = set()
    total_error = 0.0
    label = _SAMPLE_LABEL[sample]
    for i in range(n):
        j = assignment[i + 1] - 1
        if 0 <= j < m and abs(preds[i] - obs[j]) <= tolerance:
            matched_pred.add(i)
            matched_obs.add(j)
            error = abs(preds[i] - obs[j])
            total_error += error
            matches.append(BandMatchOut(
                predicted=preds[i],
                observed=obs[j],
                error=round(error, 4),
                within_tolerance=True,
                description=(f"{label}预测条带 {preds[i]:g} bp 与观测条带 "
                             f"{obs[j]:g} bp 配对，长度差 {error:g} bp，"
                             f"在容差 {tolerance:g} bp 内"),
            ))

    matches.sort(key=lambda x: x.predicted)
    missing = [preds[i] for i in range(n) if i not in matched_pred]
    extra = [obs[j] for j in range(m) if j not in matched_obs]

    return SampleIdentificationOut(
        sample=sample,
        matched_count=len(matches),
        missing_band_count=len(missing),
        extra_band_count=len(extra),
        total_error=round(total_error, 4),
        matches=matches,
        missing_bands=missing,
        extra_bands=extra,
    )


def identify(req: IdentifyRequest) -> IdentifyResponse:
    """分别对参考型、变异型做匹配，按未匹配数→总误差的顺序给出鉴别结论。"""
    reference = match_sample(
        req.reference_bands, req.observed_bands, req.tolerance, "reference",
    )
    variant = match_sample(
        req.variant_bands, req.observed_bands, req.tolerance, "variant",
    )

    ref_raw = sum(mt.error for mt in reference.matches)
    var_raw = sum(mt.error for mt in variant.matches)
    ref_unmatched = reference.missing_band_count + reference.extra_band_count
    var_unmatched = variant.missing_band_count + variant.extra_band_count
    if (ref_unmatched, ref_raw) < (var_unmatched, var_raw):
        decision = "reference"
    elif (var_unmatched, var_raw) < (ref_unmatched, ref_raw):
        decision = "variant"
    else:
        decision = "ambiguous"

    if decision == "ambiguous":
        reasoning = (
            f"参考型未匹配条带 {ref_unmatched} 条（缺失 {reference.missing_band_count}、"
            f"额外 {reference.extra_band_count}）、总误差 {reference.total_error:g} bp，"
            f"变异型未匹配条带 {var_unmatched} 条（缺失 {variant.missing_band_count}、"
            f"额外 {variant.extra_band_count}）、总误差 {variant.total_error:g} bp，"
            "两者未匹配条带数与总误差均相同，无法判定（ambiguous）"
        )
    else:
        winner, loser = (
            (reference, variant) if decision == "reference" else (variant, reference)
        )
        winner_label = _SAMPLE_LABEL[winner.sample]
        loser_label = _SAMPLE_LABEL[loser.sample]
        winner_unmatched = winner.missing_band_count + winner.extra_band_count
        loser_unmatched = loser.missing_band_count + loser.extra_band_count
        reasoning = (
            f"{winner_label}未匹配条带 {winner_unmatched} 条（缺失 "
            f"{winner.missing_band_count}、额外 {winner.extra_band_count}）、总误差 "
            f"{winner.total_error:g} bp，{loser_label}未匹配条带 {loser_unmatched} 条"
            f"（缺失 {loser.missing_band_count}、额外 {loser.extra_band_count}）、"
            f"总误差 {loser.total_error:g} bp；"
            + (
                f"{winner_label}容差内解释的条带更多，"
                if winner_unmatched != loser_unmatched
                else f"双方未匹配条带数相同，但{winner_label}总误差更小，"
            )
            + f"观测更符合{winner_label}，判定为 {winner.sample}"
        )

    return IdentifyResponse(
        tolerance=req.tolerance,
        decision=decision,
        reference=reference,
        variant=variant,
        reasoning=reasoning,
    )
