# 环状 DNA 双酶切鉴别服务

基于 FastAPI 的圆形（环状）DNA 双酶切方案推荐服务：客户端提交质粒长度、参考型与变异型的候选酶切位点、最小可见片段长度和条带分辨容差，服务枚举所有不同酶对，按圆环相邻切点计算两种样本的片段长度，将容差内的片段聚为同一条电泳条带，并推荐能够区分两种样本的酶对。

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

## 项目结构

```
.
├── app/
│   ├── main.py          # FastAPI 入口：健康检查、推荐路由、统一 422 处理
│   ├── models.py        # 请求/响应 Pydantic 模型
│   ├── validation.py    # 越界、重复位点等语义校验
│   ├── digest.py        # 圆环片段、条带聚类、酶对枚举与排序（纯函数）
│   └── config.py        # 环境变量配置（均有默认值）
├── tests/               # pytest 单元 + API 端到端测试
├── examples/request.json
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```

## 运行测试

```bash
pip install -r requirements.txt pytest httpx
python3 -m pytest -q
```
