# 环状 DNA 双酶切鉴别服务

基于 FastAPI 的圆形（环状）DNA 双酶切方案推荐与凝胶条带鉴别服务：

- **方案推荐**：客户端提交质粒长度、参考型与变异型的候选酶切位点、最小可见片段长度和条带分辨容差，服务枚举所有不同酶对，按圆环相邻切点计算两种样本的片段长度，将容差内的片段聚为同一条电泳条带，并推荐能够区分两种样本的酶对。
- **凝胶鉴别**（`POST /api/v1/digest/identify`）：客户端提交同一质粒的参考型、变异型**预测条带长度**及**一次凝胶观测长度列表**，服务在给定容差内做确定性一对一最小代价匹配，分别计算两种样本的缺失、额外条带和总误差，返回 `reference` / `variant` / `ambiguous` 判定及逐条匹配依据。

**核心功能为纯本地数学计算，不依赖数据库、缓存或任何外部服务。**

---

## 单条命令启动

需要 Docker（含 Compose 插件，Docker 20.10+ / Docker Desktop 自带）：

```bash
docker compose up --build -d
```

服务启动后监听 **`http://localhost:8000`**，两个接口（`/api/v1/digest/recommend` 与 `/api/v1/digest/identify`）均由同一服务提供，无需额外容器或环境变量。

停止：

```bash
docker compose down
```

修改代码后重建镜像：

```bash
docker compose up --build -d
```

### 不使用 Docker（本地 Python 3.11+）

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

---

## API 地址与端口

| 用途 | 方法 | 地址 |
| --- | --- | --- |
| 健康检查 | GET | `http://localhost:8000/health` |
| 双酶切方案推荐 | POST | `http://localhost:8000/api/v1/digest/recommend` |
| 凝胶观测条带鉴别 | POST | `http://localhost:8000/api/v1/digest/identify` |
| 交互式 API 文档（Swagger） | GET | `http://localhost:8000/docs` |
| OpenAPI Schema | GET | `http://localhost:8000/openapi.json` |

健康检查示例：

```bash
curl http://localhost:8000/health
# {"status":"ok","service":"circular-dna-digest","version":"1.0.0"}
```

---

## 必要配置

所有配置均有默认值，**开箱即用，无需任何外部服务**。可通过环境变量覆盖：

| 环境变量 | 默认值 | 说明 | Compose 中如何改 |
| --- | --- | --- | --- |
| `APP_HOST` | `0.0.0.0` | 服务监听地址 | 编辑 `docker-compose.yml` 的 `environment` |
| `APP_PORT` | `8000` | 容器内监听端口 | 编辑 `docker-compose.yml` 的 `environment` |
| `HOST_PORT` | `8000` | 宿主机映射端口 | `HOST_PORT=9000 docker compose up --build -d` |
| `MAX_ENZYMES` | `64` | 单次请求允许的最大酶数量（酶对数为 C(n,2)） | 编辑 `environment` |

---

## 请求与响应

### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `plasmid_length` | int | > 0，环状质粒总长度（bp），切点合法范围为 `0 ~ plasmid_length-1` |
| `min_visible_fragment` | int | ≥ 0 且 < 质粒长度；短于此长度的片段不参与条带聚类 |
| `band_resolution_tolerance` | int | ≥ 0，条带分辨容差（bp） |
| `enzymes` | array | 至少 2 种酶，酶名称不可重复 |
| `enzymes[].name` | string | 酶名称，如 `EcoRI` |
| `enzymes[].reference_sites` | int[] | 参考型候选切点（0-based，bp） |
| `enzymes[].variant_sites` | int[] | 变异型候选切点（0-based，bp） |

请求示例见 [`examples/request.json`](examples/request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/recommend \
  -H 'Content-Type: application/json' \
  -d @examples/request.json
```

### 算法说明

1. **枚举酶对**：对所有不同酶的两两组合（酶名按字典序规范化）。
2. **合并切点**：每种样本合并两酶切点；两种酶切在同一位置时合并为一个切点并记录来源酶（不视为重复位点）。
3. **圆环片段**：切点升序排列为 p0 < p1 < … < p_{n-1}，相邻切点间形成片段；末段 `(p_{n-1}, p0)` 跨过坐标原点，长度为 `p0 + L - p_{n-1}`，响应中以 `crosses_origin: true` 标记。
4. **可见性过滤**：长度 < `min_visible_fragment` 的片段不参与成带（计数为 `hidden_fragment_count`）。
5. **条带聚类**：两样本可见片段合并并按长度排序，相邻片段长度差 ≤ `band_resolution_tolerance` 归入同一条带（单连接聚类），超过则断开成新带。
6. **可区分条带**：条带内参考型与变异型的片段数量不同（`reference_only` / `variant_only` / `count_mismatch`）。
7. **推荐排序**：可区分条带数降序 → 总条带数升序 → 酶对名称字典序升序。

### 响应要点

- `recommendations`：推荐酶对，包含
  - `samples.reference` / `samples.variant`：合并后的**切点**（含来源酶）、全部**片段**及跨原点标记、可见/不可见片段计数；
  - `bands`：聚类后的**条带**（代表长度、长度范围、两样本各自片段数、成员片段、是否可区分）；
  - `differences`：**差异依据**（中文说明 + 机器可读类型）。
- `rejected_pairs`：**任一样本有效切点不足 2 个**的酶对（含原因与各自切点数量）。
- `non_discriminating_pairs`：切点足够但两样本条带完全一致、无法鉴别的酶对。
- `total_enzyme_pairs` / `recommended_count` / `rejected_pair_count` / `non_discriminating_pair_count`：汇总计数。

### 拒绝条件（HTTP 422）

以下输入问题整体拒绝，返回统一结构 `{"code":"validation_failed","errors":[...]}`，一次性报告全部错误：

- 切点越界（`< 0` 或 `≥ plasmid_length`）；
- 同一样本、同一酶的候选切点重复；
- 酶名称重复；酶少于 2 种；
- 质粒长度 ≤ 0、最小可见片段 ≥ 质粒长度、容差为负等字段约束错误。

注意：切点不足两个属于**酶对级别**的拒绝（记录在 `rejected_pairs`，HTTP 200），而非请求级别错误。

---

## 凝胶条带鉴别：`POST /api/v1/digest/identify`

对同一质粒的两种预测（参考型、变异型）与**同一次凝胶实验观测**做比对，判定观测更接近哪种样本。

### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `reference_bands` | number[] | 参考型预测条带长度（bp），允许为空列表；每项 ≥ 0 |
| `variant_bands` | number[] | 变异型预测条带长度（bp），允许为空列表；每项 ≥ 0 |
| `observed_bands` | number[] | 一次凝胶观测到的条带长度（bp），**不可为空**；每项 ≥ 0 |
| `tolerance` | number | 匹配容差（bp），必须为有限、非负数；长度差 `≤ tolerance` 才允许配对 |

请求示例见 [`examples/identify.json`](examples/identify.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/identify \
  -H 'Content-Type: application/json' \
  -d @examples/identify.json
```

### 匹配与判定算法

对参考型、变异型分别执行同一次匹配（两个样本共用 `observed_bands`，互不影响）：

1. **可行边**：预测条带 `p` 与观测条带 `o` 仅在 `|p - o| ≤ tolerance` 时可以配对；
2. **一对一约束**：每条预测条带、每条观测条带至多参与一次配对；
3. **最小代价匹配**：以绝对长度差为边权，先最大化容差内配对数，再在配对数相同的方案中最小化总绝对误差（匈牙利算法求最小权二分图匹配，复杂度 O(k³)）；
4. **缺失条带**：预测到但未匹配上观测的条带（`missing_bands`，凝胶上应有而未见）；**额外条带**：观测到但无预测可匹配的条带（`extra_bands`，凝胶上多出）；
5. **确定性**：所有长度列表先排序后匹配，结果只取决于输入数值，与列表原始顺序无关。

两种样本完成匹配后按以下优先级比较，给出 `decision`：

- 先比 **未匹配条带总数 = 缺失数 + 额外数**（容差内解释的条带越多越优）；
- 未匹配数相同，再比 **总误差 `total_error`**（全部配对绝对长度差之和）；
- 两者仍相同为 `ambiguous`。

> 说明：不能只比较配对误差——否则一个一条都没配上的样本也会以 0 误差与完全吻合的样本打平，产生退化判定。

### 响应要点

- `decision`：`reference` / `variant` / `ambiguous`；
- `reference` / `variant`：两种样本各自的匹配结果
  - `matches`：**逐条配对依据**（预测长度、观测长度、绝对误差、是否在容差内、中文说明）；
  - `missing_bands` / `extra_bands` 及对应计数；
  - `matched_count` / `total_error`；
- `reasoning`：中文判定依据（双方未匹配数与总误差的对比）。

示例请求（容差 20 bp，观测 `[102, 298, 505]`）的判定结果为 `reference`：参考型 3 条全部配对、总误差 9 bp；变异型仅 `120 ↔ 102`（误差 18）可配，缺失 `[320, 600]`、额外 `[298, 505]`。

### 拒绝条件（HTTP 422）

与推荐接口共用统一结构 `{"code":"validation_failed","errors":[...]}`：

- `observed_bands` 为空列表（预测列表允许为空，观测列表不允许）；
- 任一条带长度为负数（含预测列表与观测列表）；
- `tolerance` 为负数，或为 `NaN` / `Infinity` 等非有限值；条带长度为非有限值同样拒绝；
- 缺少必填字段、字段类型错误等 Pydantic 校验错误。

---

## 项目结构

```
.
├── app/
│   ├── main.py          # FastAPI 入口：健康检查、推荐/鉴别路由、统一 422 处理
│   ├── models.py        # 请求/响应 Pydantic 模型
│   ├── validation.py    # 越界、重复位点、非有限值等语义校验
│   ├── digest.py        # 圆环片段、条带聚类、酶对枚举与排序（纯函数）
│   ├── identify.py      # 预测/观测条带的容差内一对一最小代价匹配（纯函数）
│   └── config.py        # 环境变量配置（均有默认值）
├── tests/               # pytest 单元 + API 端到端测试
├── examples/
│   ├── request.json     # 方案推荐示例请求
│   └── identify.json    # 凝胶鉴别示例请求
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```

## 运行测试

```bash
pip install -r requirements.txt pytest httpx
python3 -m pytest -q
```
