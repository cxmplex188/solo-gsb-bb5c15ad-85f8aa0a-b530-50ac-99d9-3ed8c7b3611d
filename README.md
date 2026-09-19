# 环状 DNA 双酶切鉴别服务

基于 FastAPI 的圆形（环状）DNA 双酶切鉴别服务，提供两个互补接口：

- **方案推荐**：客户端提交质粒长度、参考型与变异型的候选酶切位点、最小可见片段长度和条带分辨容差，服务枚举所有不同酶对，按圆环相邻切点计算两种样本的片段长度，将容差内的片段聚为同一条电泳条带，并推荐能够区分两种样本的酶对。
- **样本鉴别**：客户端提交同一质粒参考型、变异型的预测条带长度和一次凝胶观测长度列表，服务在给定容差内做确定性的一对一最小代价匹配，分别计算两种样本的缺失条带、额外条带与总误差，判定观测样本为参考型、变异型或无法区分。

**核心功能为纯本地数学计算，不依赖数据库、缓存或任何外部服务。**

---

## 单条命令启动

需要 Docker（含 Compose 插件，Docker 20.10+ / Docker Desktop 自带）：

```bash
docker compose up --build -d
```

服务启动后监听 **`http://localhost:8000`**。

停止：

```bash
docker compose down
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
| 凝胶条带样本鉴别 | POST | `http://localhost:8000/api/v1/digest/identify` |
| 交互式 API 文档（Swagger） | GET | `http://localhost:8000/docs` |
| OpenAPI Schema | GET | `http://localhost:8000/openapi.json` |

健康检查示例：

```bash
curl http://localhost:8000/health
# {"status":"ok","service":"circular-dna-digest","version":"1.1.0"}
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
| `MAX_BANDS` | `512` | identify 接口每个条带列表允许的最大条带数 | 编辑 `environment` |

---

## 请求与响应

### 双酶切方案推荐（POST `/api/v1/digest/recommend`）

#### 请求体（JSON）

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

#### 算法说明（推荐）

1. **枚举酶对**：对所有不同酶的两两组合（酶名按字典序规范化）。
2. **合并切点**：每种样本合并两酶切点；两种酶切在同一位置时合并为一个切点并记录来源酶（不视为重复位点）。
3. **圆环片段**：切点升序排列为 p0 < p1 < … < p_{n-1}，相邻切点间形成片段；末段 `(p_{n-1}, p0)` 跨过坐标原点，长度为 `p0 + L - p_{n-1}`，响应中以 `crosses_origin: true` 标记。
4. **可见性过滤**：长度 < `min_visible_fragment` 的片段不参与成带（计数为 `hidden_fragment_count`）。
5. **条带聚类**：两样本可见片段合并并按长度排序，相邻片段长度差 ≤ `band_resolution_tolerance` 归入同一条带（单连接聚类），超过则断开成新带。
6. **可区分条带**：条带内参考型与变异型的片段数量不同（`reference_only` / `variant_only` / `count_mismatch`）。
7. **推荐排序**：可区分条带数降序 → 总条带数升序 → 酶对名称字典序升序。

#### 响应要点（推荐）

- `recommendations`：推荐酶对，包含
  - `samples.reference` / `samples.variant`：合并后的**切点**（含来源酶）、全部**片段**及跨原点标记、可见/不可见片段计数；
  - `bands`：聚类后的**条带**（代表长度、长度范围、两样本各自片段数、成员片段、是否可区分）；
  - `differences`：**差异依据**（中文说明 + 机器可读类型）。
- `rejected_pairs`：**任一样本有效切点不足 2 个**的酶对（含原因与各自切点数量）。
- `non_discriminating_pairs`：切点足够但两样本条带完全一致、无法鉴别的酶对。
- `total_enzyme_pairs` / `recommended_count` / `rejected_pair_count` / `non_discriminating_pair_count`：汇总计数。

### 凝胶条带样本鉴别（POST `/api/v1/digest/identify`）

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `reference_bands` | number[] | 参考型预测条带长度（bp），每条 ≥ 0，可为空列表 |
| `variant_bands` | number[] | 变异型预测条带长度（bp），每条 ≥ 0，可为空列表 |
| `observed_bands` | number[] | 一次凝胶观测到的条带长度（bp），每条 ≥ 0，**至少 1 条** |
| `band_match_tolerance` | number | ≥ 0，条带匹配容差（bp） |

请求示例见 [`examples/identify-request.json`](examples/identify-request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/identify \
  -H 'Content-Type: application/json' \
  -d @examples/identify-request.json
```

#### 匹配与判定规则

1. **一对一最小代价匹配**：对每种样本，将其预测条带与观测条带分别按长度升序排列后做一对一配对；仅当 `|预测长度 − 观测长度| ≤ band_match_tolerance` 时允许配对。目标为先**最大化匹配对数**，再**最小化总误差**（各配对长度差之和）。一维绝对差代价满足 Monge 性质，最优解必可取不交叉匹配，服务用序列比对式动态规划在 O(n·m) 内求得；并列时按固定优先级选择（等效于较短条带优先配对），结果完全确定且与列表提交顺序无关。
2. **缺失 / 额外条带**：未匹配到观测条带的预测条带为**缺失条带**（`missing_bands`），未匹配到任何预测条带的观测条带为**额外条带**（`extra_bands`）。
3. **判定**：比较两种样本的未匹配条带数（缺失 + 额外），少者判为对应样本（`reference` / `variant`）；并列时总误差小者胜出；两者完全一致时判定为 `ambiguous`（无法区分）。

#### 响应要点（鉴别）

- `samples.reference` / `samples.variant`：每种样本的
  - `matched_count` / `missing_count` / `extra_count` / `total_error`：匹配对数、缺失数、额外数与总误差（bp，保留两位小数）；
  - `matches`：**逐条匹配依据**（预测长度、观测长度、长度差及中文说明），按预测长度升序；
  - `missing_bands` / `extra_bands`：缺失、额外条带长度列表（升序）。
- `verdict`：`reference` / `variant` / `ambiguous`；`verdict_reason` 为中文判定依据，可完全由响应中的计数与总误差复算。

### 拒绝条件（HTTP 422）

以下输入问题整体拒绝，返回统一结构 `{"code":"validation_failed","errors":[...]}`，一次性报告全部错误：

推荐接口（`/api/v1/digest/recommend`）：

- 切点越界（`< 0` 或 `≥ plasmid_length`）；
- 同一样本、同一酶的候选切点重复；
- 酶名称重复；酶少于 2 种；
- 质粒长度 ≤ 0、最小可见片段 ≥ 质粒长度、容差为负等字段约束错误。

鉴别接口（`/api/v1/digest/identify`）：

- 观测条带列表为空（`observed_bands` 至少 1 条）；
- 任一条带长度（参考 / 变异 / 观测）为负数，或匹配容差为负、缺失、非数值；
- 任一条带列表超过 `MAX_BANDS` 上限。

注意：推荐接口中切点不足两个属于**酶对级别**的拒绝（记录在 `rejected_pairs`，HTTP 200），而非请求级别错误。

---

## 项目结构

```
.
├── app/
│   ├── main.py          # FastAPI 入口：健康检查、推荐/鉴别路由、统一 422 处理
│   ├── models.py        # 请求/响应 Pydantic 模型
│   ├── validation.py    # 越界、重复位点、条带规模等语义校验
│   ├── digest.py        # 圆环片段、条带聚类、酶对枚举与排序（纯函数）
│   ├── identify.py      # 条带一对一最小代价匹配与样本判定（纯函数）
│   └── config.py        # 环境变量配置（均有默认值）
├── tests/               # pytest 单元 + API 端到端测试
├── examples/request.json
├── examples/identify-request.json
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```

## 运行测试

```bash
pip install -r requirements.txt pytest httpx
python3 -m pytest -q
```
