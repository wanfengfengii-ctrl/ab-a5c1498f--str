# 混合 DNA 样本裁决系统

法医检验室复核混合 DNA 样本时，需要从候选供样人的 STR 等位基因型中还原**最少供样人组合及其整数模板份数**。本系统在「候选人子集 × 正整数份数」的**联合空间**中裁决，避免按单个位点各自猜测而产生无法同时解释全样本的结论。

## 裁决模型（可复核约定）

- 录入 4–8 名唯一编号候选人、3–6 个位点，每人每位点两条等位基因；每个位点按等位基因汇总观测拷贝数，并给定统一的非负整数误差上限。
- 设各位点共同的观测总深度为 `D`（各位点观测拷贝数之和必须相等，否则视为非法输入）。
- 从候选人中选取 1–3 人，并分配**正整数**份数 `n_i`，约束 `Σ n_i = D`。
- 每份模板在该位点的两条等位基因上各贡献 1/2 份（纯合子合计 1 份），即等位基因 `a` 的预测拷贝数
  `pred(a) = Σ_i n_i · m_i(a) / 2`，其中 `m_i(a) ∈ {0, 1, 2}` 为候选人 `i` 携带 `a` 的条数。
- **可行性**：每个位点每种等位基因（观测与候选基因型的并集）都满足 `|pred(a) − obs(a)| ≤ 误差上限`。
- **目标**：先最少化供样人数，再最小化全部位点全部等位基因的绝对残差和 `Σ|pred − obs|`。
- **唯一性**：同一最优值下，若按候选人编号排序的份数向量只有一种，判定为 `unique`；否则判定为 `ambiguous` 并列出两份不同见证。无可行组合时为 `no_solution`。非法输入返回 422。歧义、无解、非法输入均不会伪装为唯一解。

## 快速开始（Docker）

```bash
# 启动服务（默认宿主机端口 8080，可用 HOST_PORT 覆盖）
docker compose up --build app
HOST_PORT=9000 docker compose up --build app

# 打开前端
open http://localhost:8080/

# 一次性验证：构建镜像 → 运行代码测试 → 对运行中的服务做 HTTP 冒烟 → 自动退出并以退出码报告结果
docker compose up --build --abort-on-container-exit --exit-code-from verify verify
echo $?   # 0 = 全部通过，非 0 = 存在失败项
docker compose down
```

- `app` 服务自带健康检查（`GET /api/health`，Dockerfile 与 Compose 双重定义）。
- `verify` 为一次性服务：等待 `app` 健康后依次执行 `pytest` 与 `verify/smoke.py`，随后自行退出。

## 本地开发

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/uvicorn app.main:app --reload --port 8000

# 测试与冒烟
.venv/bin/python -m pytest tests -q
APP_URL=http://127.0.0.1:8000 .venv/bin/python verify/smoke.py
```

## API

### `POST /api/mixture-adjudications`

请求体：

```json
{
  "tolerance": 0,
  "candidates": [
    {"id": "C1", "genotypes": {"L1": ["a", "b"], "L2": ["x", "y"], "L3": ["p", "q"]}},
    {"id": "C2", "genotypes": {"L1": ["a", "c"], "L2": ["x", "x"], "L3": ["p", "p"]}},
    {"id": "C3", "genotypes": {"L1": ["b", "c"], "L2": ["y", "z"], "L3": ["q", "r"]}},
    {"id": "C4", "genotypes": {"L1": ["c", "c"], "L2": ["z", "z"], "L3": ["r", "r"]}}
  ],
  "loci": [
    {"name": "L1", "observed": {"a": 1, "b": 2, "c": 1}},
    {"name": "L2", "observed": {"x": 1, "y": 2, "z": 1}},
    {"name": "L3", "observed": {"p": 1, "q": 2, "r": 1}}
  ]
}
```

响应（唯一解示例，`predicted` 可能为半份）：

```json
{
  "status": "unique",
  "totalDepth": 4,
  "tolerance": 0,
  "objective": {"contributors": 2, "absoluteResidual": 0.0},
  "solution": {
    "contributors": [{"id": "C1", "copies": 2}, {"id": "C3", "copies": 2}],
    "absoluteResidual": 0.0,
    "perLocus": [
      {"locus": "L1", "residual": 0.0,
       "alleles": [{"allele": "a", "observed": 1, "predicted": 1.0, "absDiff": 0.0}]}
    ]
  },
  "witnesses": null,
  "examined": {"subsets": 10, "assignments": 22}
}
```

- `status`：`unique` / `ambiguous` / `no_solution`。
- `ambiguous` 时 `witnesses` 含两份不同的最优份数向量；`unique` 时 `solution` 为唯一最优解；`no_solution` 时两者均为 `null`。
- 非法输入（人数/位点数越界、编号重复、总深度不一致、负拷贝数等）返回 `422` 及明细。

### `GET /api/health`

返回 `{"status": "ok"}`，用于容器健康检查。

## 工程约束

- 观测总深度上限 500、每位点等位基因并集上限 40（保证组合枚举在有限时间内完成；最坏规模约 5 秒）。
- 求解器按二倍拷贝数做整数运算，避免半份带来的浮点误差。

## 目录结构

```
app/
  main.py            # FastAPI 入口：/api/health、/api/mixture-adjudications、前端页面
  schemas.py         # 请求结构校验（非法输入 → 422）
  solver.py          # 联合空间裁决求解器（子集 × 正整数份数，numpy 向量化枚举）
  static/index.html  # 录入与裁决结果展示页面（原生 JS）
tests/               # pytest：求解器单元测试 + API 测试
verify/
  run.sh             # 一次性验证入口：pytest → HTTP 冒烟
  smoke.py           # HTTP 冒烟脚本（健康检查、页面、唯一/歧义/422）
Dockerfile           # 单镜像：运行服务 + 测试/冒烟，含 HEALTHCHECK
docker-compose.yml   # app（可配置 HOST_PORT）+ verify（一次性）
```
