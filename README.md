# 混合 DNA 样本供样人裁决（STR Mixture Adjudicator）

法医检验室复核混合 DNA 样本时使用的全栈小工具：检验员在页面录入候选供样人的
STR 等位基因型与混合样本的观测拷贝数，后端在**完整候选人子集 × 正整数份数**的
联合空间中裁决，前端展示可复核的逐位点明细。非法输入、无解与歧义均不会被
伪装成唯一解。

## 裁决模型

- 录入 **4–8 名**唯一编号候选人、**3–6 个**位点；每名候选人在每个位点有恰好
  两条等位基因（可相同，即纯合子）。
- 每个位点按等位基因汇总观测拷贝数；所有位点的**观测总深度必须一致**，记为 D。
- 选中供样人数限 **1–3 人**，每人一个正整数模板份数，且**份数之和等于 D**。
- 供样人 i（份数 cᵢ）对位点某等位基因 a 的预测贡献为
  `cᵢ × (a 在其基因型中的剂量 0/1/2) / 2`，因此每个位点的预测总深度恰为 D。
- 对每个位点、每种相关等位基因（观测中出现或选中者基因型中出现），要求
  `|预测 − 观测| ≤ 统一误差上限`。
- 目标：先**最少化供样人数**，再**最小化全部绝对残差之和**。
- 唯一性：达到同一最优值的“按候选人编号排序的份数向量”若只有一种，判定
  `unique`；否则判定 `ambiguous` 并列出**两份不同见证**；不存在可行组合时判定
  `no_solution`。

## 目录结构

```
app/
  main.py          FastAPI 入口（静态页面 + API）
  models.py        请求体校验（非法输入一律 422）
  solver.py        联合空间穷举求解器（numpy 向量化）
  static/          前端页面（index.html / app.js / style.css）
tests/             pytest 单元测试与 API 集成测试
scripts/smoke.py   HTTP 冒烟脚本（verify 服务使用）
Dockerfile         单镜像：运行服务 + 内置健康检查
docker-compose.yml web 服务（可配置宿主机端口）+ verify 一次性服务
```

## 本地运行

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
# 打开 http://localhost:8000 ，可点击“载入示例”后“裁决”
```

运行测试：

```bash
python -m pytest            # 代码测试
python scripts/smoke.py     # 对运行中的服务做 HTTP 冒烟（BASE_URL 可覆盖）
```

## Docker 运行

```bash
docker compose up --build            # 默认宿主机端口 8000
HOST_PORT=9000 docker compose up --build   # 自定义宿主机端口
```

- `web` 服务：Dockerfile 与 Compose 均配置了 `/api/health` 健康检查。
- 宿主机端口通过环境变量 `HOST_PORT` 配置（默认 8000），容器内端口固定 8000。

## 一键验证（构建 + 代码测试 + HTTP 冒烟）

```bash
docker compose up --build --exit-code-from verify --abort-on-container-exit verify
echo $?   # 0 = 全部通过；非 0 = 失败
```

`verify` 是一次性服务：镜像构建完成后，它等待 `web` 健康，先运行
`pytest` 代码测试，再对 `web` 执行 `scripts/smoke.py` 的 HTTP 冒烟
（健康检查、页面、唯一解 / 歧义 / 无解 / 非法输入四类场景），随后自行退出，
并以退出码报告结果。

## API

### `POST /api/mixture-adjudications`

请求体：

```json
{
  "candidates": [
    {"id": "C1", "genotypes": {"L1": ["12", "13"], "L2": ["7", "9"], "L3": ["16", "17"]}},
    {"id": "C2", "genotypes": {"L1": ["11", "12"], "L2": ["8", "8"], "L3": ["15", "18"]}},
    {"id": "C3", "genotypes": {"L1": ["14", "14"], "L2": ["9", "10"], "L3": ["17", "18"]}},
    {"id": "C4", "genotypes": {"L1": ["13", "15"], "L2": ["7", "10"], "L3": ["16", "16"]}},
    {"id": "C5", "genotypes": {"L1": ["11", "15"], "L2": ["8", "11"], "L3": ["15", "15"]}}
  ],
  "loci": ["L1", "L2", "L3"],
  "observed": {
    "L1": {"12": 2, "13": 2, "14": 2},
    "L2": {"7": 2, "9": 3, "10": 1},
    "L3": {"16": 2, "17": 3, "18": 1}
  },
  "tolerance": 0
}
```

响应（`status` 为 `unique` / `ambiguous` / `no_solution`；非法输入返回 422）：

```json
{
  "status": "unique",
  "totalDepth": 6,
  "tolerance": 0,
  "objective": {"contributors": 2, "residualSum": 0},
  "optimalVectorCount": 1,
  "solution": {
    "contributors": [{"candidateId": "C1", "copies": 4}, {"candidateId": "C3", "copies": 2}],
    "residualSum": 0,
    "loci": [
      {"locus": "L1", "alleles": [
        {"allele": "12", "observed": 2, "predicted": 2, "residual": 0}
      ]}
    ]
  },
  "witnesses": null
}
```

- `ambiguous`：`solution` 为 `null`，`witnesses` 给出两份不同见证，
  `optimalVectorCount` 报告并列最优的份数向量总数。
- `no_solution`：`objective` 为 `null`，`witnesses` 为空数组。
- 预测值与残差可能为 `x.5`（奇数份数的杂合子），JSON 数值精确表示。

### `GET /api/health`

返回 `{"status": "ok"}`，供 Docker 健康检查使用。

## 规模限制

- 候选人 4–8 名、位点 3–6 个、每点位观测等位基因 ≤ 32 种。
- 观测总深度 1–1000（份数组合穷举随深度平方增长，深度 1000 时最坏情形约数秒）。
- 误差上限 0–1000。
