# 化学反应计量平衡服务

基于 **Python + FastAPI** 的化学反应方程式配平微服务。客户端提交包含
**元素、括号与电荷**的反应物 / 生成物分子式，服务端：

1. 用递归下降解析器解析化学式（支持三种括号、嵌套、下标乘数、离子电荷、电子 `e-`）；
2. 建立 **元素守恒 + 总电荷守恒** 矩阵（反应物列取正、生成物列取负）；
3. 用 `fractions.Fraction` 做精确 **高斯-若尔当消元（RREF）** 求零空间，
   把唯一的有理数解放大、约分得到 **最小互质正整数系数**；
4. `/balance/variants` 在自由变量不超过 3 个时，枚举零空间内所有
   `1..max_coefficient` 的互质正整数配平方案；
5. `/stoichiometry/batch` 在唯一配平之上，用 `Decimal` 精确计算反应进度、
   并列限量试剂、各反应物消耗/剩余与各生成物理论产量；
6. `/stoichiometry/plan` 做反向投料规划：按一个或多个生成物的目标摩尔数，
   取各目标量/配平系数的最大值得到满足全部目标的最小反应进度，给出各反应物
   理论需求量、各目标产物达成量与超额量；
7. 对同一规范化反应做幂等缓存，重复请求返回相同的 `reaction_id` 与结果。

---

## 1. 一条命令启动

需要已安装 Docker 与 Docker Compose v2：

```bash
docker compose up --build -d
```

服务监听 **宿主机端口 `8000`**（映射容器 8000）。

| 入口 | 地址 |
| --- | --- |
| 健康检查 | <http://localhost:8000/health> |
| Swagger 交互文档 | <http://localhost:8000/docs> |
| ReDoc 文档 | <http://localhost:8000/redoc> |

> 想换宿主机端口：编辑 `docker-compose.yml` 的 `ports`，例如
> `"9000:8000"` 后访问 9000 端口。容器内固定监听 8000。

停止服务：

```bash
docker compose down
```

启动后可直接调用批量计量接口（宿主机端口 **8000**）：

```bash
curl -X POST http://localhost:8000/api/v1/stoichiometry/batch \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["CH4","O2"],"products":["CO2","H2O"],
       "amounts":[{"species":"CH4","moles":"1"},{"species":"O2","moles":"3"}]}'
```

反向投料规划——要求至少产出 4 mol NH₃，看需要多少 N₂ / H₂：

```bash
curl -X POST http://localhost:8000/api/v1/stoichiometry/plan \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["N2","H2"],"products":["NH3"],
       "targets":[{"species":"NH3","moles":"4"}]}'
```

### 不用 Docker 的本地运行方式

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

运行测试：

```bash
pip install -r requirements-dev.txt
pytest          # 175 个测试
```

---

## 2. API 一览（前缀 `/api/v1`）

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/api/v1/balance` | 配平反应，返回最小正整数系数与电荷守恒结果 |
| POST | `/api/v1/balance/variants` | 枚举自由变量 ≤3 时所有互质正整数配平方案 |
| POST | `/api/v1/analyze` | 逐元素差异说明：配平前（单位系数）与配平后两个视图 |
| POST | `/api/v1/stoichiometry/batch` | 投料摩尔数批量计量：反应进度、并列限量试剂、消耗/剩余与理论产量 |
| POST | `/api/v1/stoichiometry/plan` | 反向投料规划：按生成物目标摩尔数求最小反应进度、理论投料需求与达成/超额量 |
| GET | `/api/v1/reactions/{reaction_id}` | 按幂等 ID 取回之前的结果 |
| GET | `/health` | 健康检查 |

### 2.1 配平 `POST /api/v1/balance`

请求体：

```json
{
  "reactants": ["KMnO4", "HCl"],
  "products": ["KCl", "MnCl2", "Cl2", "H2O"]
}
```

响应 `200`：

```json
{
  "reaction_id": "rxn_7ed134d0c6ec5699",
  "cached": false,
  "normalized": {
    "reactants": ["HCl", "KMnO4"],
    "products": ["KCl", "Cl2", "MnCl2", "H2O"],
    "equation": "16HCl + 2KMnO4 -> 2KCl + 5Cl2 + 2MnCl2 + 8H2O"
  },
  "coefficients": {
    "reactants": [
      {"input": "HCl", "normalized": "HCl", "coefficient": 16},
      {"input": "KMnO4", "normalized": "KMnO4", "coefficient": 2}
    ],
    "products": [
      {"input": "KCl", "normalized": "KCl", "coefficient": 2},
      {"input": "Cl2", "normalized": "Cl2", "coefficient": 5},
      {"input": "MnCl2", "normalized": "MnCl2", "coefficient": 2},
      {"input": "H2O", "normalized": "H2O", "coefficient": 8}
    ],
    "vector": [16, 2, 2, 5, 2, 8]
  },
  "balanced": true,
  "elements": ["Cl", "H", "K", "Mn", "O"],
  "charge_balance": {"reactant_total": 0, "product_total": 0, "difference": 0, "balanced": true}
}
```

`curl`：

```bash
curl -X POST http://localhost:8000/api/v1/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["KMnO4","HCl"],"products":["KCl","MnCl2","Cl2","H2O"]}'
```

### 2.2 逐元素差异 `POST /api/v1/analyze`

返回两个视图：

* `raw`：所有物种按系数 1 统计（即“配平前”的差异，便于解释哪里不平衡）；
* `balanced`：用求得的系数统计，此时每个元素与电荷的 `difference` 均为 0。

以高锰酸根酸性还原半反应为例：

```bash
curl -X POST http://localhost:8000/api/v1/analyze \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["MnO4-","H+","e-"],"products":["Mn2+","H2O"]}'
```

`raw.charge` 为 `{"reactant_total": -1, "product_total": 2, "difference": -3}`，
`balanced` 视图给出 `MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O`，电荷两侧均为 +2。

即使反应**无法唯一配平**，`/analyze` 也返回 `200`，同时带 `raw` 视图与
`error` 字段（`balanced_ok=false`、`balanced=null`），方便客户端展示诊断信息。

### 2.3 枚举所有配平方案 `POST /api/v1/balance/variants`

`/balance` 只返回唯一配平；当守恒方程存在自由变量（如不完整的燃烧反应）时，
本接口在零空间中枚举**所有**满足条件的正整数方案：

* 每个系数都在 `1 ~ max_coefficient`（**2～30**，必填）；
* 元素与总电荷守恒（由零空间解天然保证）；
* 系数最大公约数为 1（只保留本原解，`2H2+O2→2H2O` 与 `4H2+2O2→4H2O`
  视为同一方案）。

零空间自由变量（= 零空间维数 − 1）超过 **3** 个（即维数 > 3）时，枚举空间
过大，返回 `422 VARIANT_SPACE_TOO_LARGE`。结果先按**系数总和**升序、再按
**规范化物种顺序的系数向量字典序**排序，取前 `limit`（**1～100**，必填）条，
并用 `truncated` 标识是否发生截取。唯一配平的反应同样适用，返回只含单个
方案的列表。

```bash
curl -X POST http://localhost:8000/api/v1/balance/variants \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["C","O2"],"products":["CO","CO2"],"max_coefficient":10,"limit":5}'
```

`200` 响应：

```json
{
  "reaction_id": "rxn_aa766b24597bded4",
  "normalized": {
    "reactants": ["C", "O2"],
    "products": ["CO", "CO2"],
    "equation": "C + O2 -> CO + CO2"
  },
  "max_coefficient": 10,
  "limit": 5,
  "nullity": 2,
  "count": 5,
  "total": 15,
  "truncated": true,
  "variants": [
    {
      "coefficients": {
        "reactants": [
          {"input": "C", "normalized": "C", "coefficient": 3},
          {"input": "O2", "normalized": "O2", "coefficient": 2}
        ],
        "products": [
          {"input": "CO", "normalized": "CO", "coefficient": 2},
          {"input": "CO2", "normalized": "CO2", "coefficient": 1}
        ],
        "vector": [3, 2, 2, 1]
      },
      "sum": 8,
      "equation": "3C + 2O2 -> 2CO + CO2"
    }
  ]
}
```

字段说明：

| 字段 | 含义 |
| --- | --- |
| `nullity` | 守恒矩阵零空间维数（本接口只返回 1～3） |
| `total` | 在 `max_coefficient` 范围内枚举出的方案总数（截取前） |
| `count` | 本次实际返回的条数（`total` 与 `limit` 的较小值） |
| `truncated` | `total > limit` 时为 `true` |
| `variants[].sum` | 该方案系数总和，即第一排序键 |
| `variants[].coefficients.vector` | 按规范化后 reactants + products 顺序的系数向量 |

不存在正整数解时（如唯一解要求某物种系数为 0）返回 `200` 且
`variants=[]`、`total=0`；守恒本身无解返回 `422 IMPOSSIBLE_BALANCE`。
`reaction_id` 与 `/balance`、`/analyze` 采用同一规范化签名（变体结果不写
入幂等缓存）。

### 2.4 批量限量试剂计量 `POST /api/v1/stoichiometry/batch`

在唯一配平的基础上，给出每种**反应物的投料摩尔数**，服务端以精确十进制
（`Decimal`，全程不经过 `float`）计算：

* 反应进度 ξ = min(**投料量 ÷ 配平系数**)；
* **全部并列限量试剂**——当投料恰好按系数比时可能多于一种（全部投料为 0 时全部限量）；
* 每种反应物的投料量 `fed`、消耗量 `consumed = ξ·系数` 与剩余量
  `remaining = fed − consumed`；
* 每种生成物的理论产量 `theoretical_yield = ξ·系数`。

请求体在 `reactants` / `products` 之外增加 `amounts`，**每种规范化反应物
恰好一条记录**（等价写法可识别，但不能缺失、重复或出现未知物种）；
摩尔数必须是**非负十进制字符串**（不要传 JSON 数字，以避免浮点误差）：

```bash
curl -X POST http://localhost:8000/api/v1/stoichiometry/batch \
  -H 'Content-Type: application/json' \
  -d '{
    "reactants": ["CH4", "O2"],
    "products": ["CO2", "H2O"],
    "amounts": [
      {"species": "CH4", "moles": "1"},
      {"species": "O2",  "moles": "3"}
    ]
  }'
```

`200` 响应（十进制值均为字符串；物种按规范化顺序排列，调换书写顺序结果不变）：

```json
{
  "reaction_id": "rxn_……",
  "cached": false,
  "normalized": {
    "reactants": ["CH4", "O2"],
    "products": ["CO2", "H2O"],
    "equation": "CH4 + 2O2 -> CO2 + 2H2O"
  },
  "reaction_extent": {"value": "1", "unit": "mol"},
  "limiting_reactants": ["CH4"],
  "reactants": [
    {"input": "CH4", "normalized": "CH4", "coefficient": 1,
     "fed": "1", "consumed": "1", "remaining": "0", "is_limiting": true},
    {"input": "O2", "normalized": "O2", "coefficient": 2,
     "fed": "3", "consumed": "2", "remaining": "1", "is_limiting": false}
  ],
  "products": [
    {"input": "CO2", "normalized": "CO2", "coefficient": 1, "theoretical_yield": "1"},
    {"input": "H2O", "normalized": "H2O", "coefficient": 2, "theoretical_yield": "2"}
  ]
}
```

投料比恰好等于系数比（如 `CH4=2.5, O2=5`）时 `limiting_reactants` 同时返回
`["CH4", "O2"]`，二者 `remaining` 都是 `"0"`。

**无限小数**（分母含 2、5 以外的因子，如 ξ = 1/7）按 **28 位有效数字**
银行家舍入输出，例如 `"0.1428571428571428571428571429"`；有限小数给出精确值。

错误沿用统一结构化格式：

| 场景 | HTTP / code |
| --- | --- |
| 反应式或投料物种无法解析 | `400 FORMULA_PARSE_ERROR` |
| 反应式同侧物种重复 | `409 DUPLICATE_SPECIES` |
| 投料缺失、重复、含未知物种（详情在 `details.missing/duplicate/unknown`） | `422 INVALID_REQUEST` |
| 摩尔数非法：负数、`NaN`、`Infinity`、非数字字符串或 JSON 数字（详情 `details.invalid_moles`） | `422 INVALID_REQUEST` |
| 反应无法唯一配平 / 无法配平 | `422 NON_UNIQUE_BALANCE` / `IMPOSSIBLE_BALANCE` 等 |

`reaction_id` 与 `/balance`、`/analyze` 共享同一规范化签名——同一反应调换
反应物 / 生成物书写顺序或使用等价化学式，得到相同的 `reaction_id` 与完全
一致的计算结果（仅 `input` 字段回显原始写法）。

### 2.5 反向投料规划 `POST /api/v1/stoichiometry/plan`

与 `/stoichiometry/batch` 方向相反：客户端为**一个或多个生成物**提交目标
摩尔数（允许只覆盖部分生成物），服务端以精确十进制（全程不经过 `float`）
计算：

* 满足全部目标所需的**最小反应进度** ξ = max(**目标量 ÷ 配平系数**)
  （取最大值是为了让每个目标产物的达成量都不低于其目标量）；
* 全部并列取得最大值的**驱动目标产物** `driving_targets`——目标比恰好
  等于系数比时可多于一种；
* 每种反应物的**理论需求量** `required = ξ·系数`；
* 每个提交目标的产物：目标量 `target`、达成量 `achieved = ξ·系数`、
  **超额量** `excess = achieved − target`（驱动产物为 `"0"`，其余为正）；
* 全部生成物（含未提交目标的）在该进度下的达成量。

请求体用 `targets` 代替 `amounts`；目标物种必须是规范化生成物之一
（等价写法可识别），不能重复或出现未知物种，摩尔数同样必须是非负十进制
字符串：

```bash
curl -X POST http://localhost:8000/api/v1/stoichiometry/plan \
  -H 'Content-Type: application/json' \
  -d '{
    "reactants": ["NH3", "O2"],
    "products": ["NO", "H2O"],
    "targets": [
      {"species": "NO",  "moles": "1"},
      {"species": "H2O", "moles": "6"}
    ]
  }'
```

`4NH3 + 5O2 -> 4NO + 6H2O` 中两个目标比值分别为 1/4 与 6/6，ξ 取最大值
`1`（由 `H2O` 驱动），因此 `NO` 实际达成 4 mol、超额 3 mol：

```json
{
  "reaction_id": "rxn_……",
  "cached": false,
  "normalized": {
    "reactants": ["NH3", "O2"],
    "products": ["NO", "H2O"],
    "equation": "4NH3 + 5O2 -> 4NO + 6H2O"
  },
  "reaction_extent": {"value": "1", "unit": "mol"},
  "driving_targets": ["H2O"],
  "reactants": [
    {"input": "NH3", "normalized": "NH3", "coefficient": 4, "required": "4"},
    {"input": "O2",  "normalized": "O2",  "coefficient": 5, "required": "5"}
  ],
  "targets": [
    {"input": "NO", "normalized": "NO", "coefficient": 4,
     "target": "1", "achieved": "4", "excess": "3", "is_driving": false},
    {"input": "H2O", "normalized": "H2O", "coefficient": 6,
     "target": "6", "achieved": "6", "excess": "0", "is_driving": true}
  ],
  "products": [
    {"input": "NO", "normalized": "NO", "coefficient": 4, "achieved": "4"},
    {"input": "H2O", "normalized": "H2O", "coefficient": 6, "achieved": "6"}
  ]
}
```

单目标示例：`3H2 + N2 -> 2NH3`，目标 `NH3 = 4` mol 时 ξ = 2，理论需求
`H2 = 6` mol、`N2 = 2` mol，达成量恰为 4 mol、超额为 0。目标摩尔数为
`"0"` 时 ξ 与所有理论需求均为 `"0"`。无限小数同样按 **28 位有效数字**
银行家舍入输出（如目标 `CO2 = 1` mol 于 `Fe2O3 + 3CO -> 2Fe + 3CO2`
时 ξ = 1/3 → `"0.3333333333333333333333333333"`）。

错误结构与 `/stoichiometry/batch` 一致：

| 场景 | HTTP / code |
| --- | --- |
| 反应式或目标物种无法解析 | `400 FORMULA_PARSE_ERROR` |
| 反应式同侧物种重复 | `409 DUPLICATE_SPECIES` |
| 目标重复或包含非生成物物种（详情 `details.duplicate/unknown`；允许只提交部分生成物） | `422 INVALID_REQUEST` |
| 摩尔数非法：负数、`NaN`、`Infinity`、非数字字符串或 JSON 数字（详情 `details.invalid_moles`） | `422 INVALID_REQUEST` |
| 反应无法唯一配平 / 无法配平 | `422 NON_UNIQUE_BALANCE` / `IMPOSSIBLE_BALANCE` 等 |

### 2.6 幂等与查询

* 幂等键与书写顺序、括号风格、空白以及电荷的等价写法无关，只取决于每一侧
  物种的 **（元素组成多重集, 电荷）**。
* 因此下面两个请求得到 **同一个 `reaction_id`**，第二次 `cached=true`：

```bash
curl -X POST .../api/v1/balance -d '{"reactants":["Ca(OH)2","HCl"],"products":["CaCl2","H2O"]}'
curl -X POST .../api/v1/balance -d '{"reactants":["Ca [OH]2"," HCl"],"products":["CaCl2","H2 O"]}'
```

* 用 `GET /api/v1/reactions/{reaction_id}` 可随时取回缓存结果（LRU，默认 1 万条）。

---

## 3. 化学式书写约定

| 特性 | 写法示例 |
| --- | --- |
| 元素（118 种） | `H`、`Fe`、`O`、`Uuo`（必须大写开头） |
| 下标 | `H2O`、`Al2(SO4)3` |
| 三种括号（等价，可嵌套） | `Ca(OH)2`、`K4[Fe(CN)6]`、`Mg3{Fe(CN)6}2` |
| 显式电荷 | `SO4^2-`、`Fe^3+`、`NH4^+`、`O^{2-}`、`[Cu(NH3)4]^2+` |
| 纯文本电荷 | `Fe3+`、`O2-`、`NH4+`、`MnO4-`、`OH-` |
| 电荷前空格消歧 | `SO4 2-` 等价于 `SO4^2-` |
| 电子 / 正电子 | `e-`、`e`（默认 −1）、`e+` |
| 空白 | 化学式内部空白一律忽略 |

**末尾 `数字±`（无 `^`）消歧规则**（沿用常见化学软件约定）：

* 物种只含 **一个原子/基团** 时，数字是电荷：`Fe3+` → Fe³⁺、`O2-` → O²⁻；
* 含 **多个原子** 时，数字是下标、符号单独表示 ±1：`NH4+` → NH₄⁺、`MnO4-` → MnO₄⁻；
* 数字前有空格时一定是电荷：`SO4 2-` → SO₄²⁻；
* 任何不确定的场景都建议写显式 `^`，无歧义。

---

## 4. 配平原理

设所有物种系数向量为 `x`（未知）。对每种元素写一行守恒方程，再加一行
总电荷守恒；反应物一侧的列取正，生成物一侧的列取负，得到齐次方程组
`A x = 0`。

对整数矩阵的每个元素套上 `Fraction` 做 RREF：

* **零空间维数 = 1**：基向量所有分量非零且同号 → 整体取正，按分母最小公倍数
  放大、再除以最大公约数，得到最小互质正整数系数；
* **维数 = 0**：`IMPOSSIBLE_BALANCE`（元素/电荷守恒本身无解）；
* **维数 > 1**：`/balance` 返回 `NON_UNIQUE_BALANCE`（自由度 = 维数 − 1，
  常见于缺少电子、H₂O、H⁺/OH⁻ 等物种）；
* 维数为 1 但有零分量 → `ZERO_COEFFICIENT`；符号不一致 → `INVALID_STOICHIOMETRY`。

**多变体枚举（`/balance/variants`）**：RREF 后第 k 个基向量在第 k 个自由列
上为 1、其余自由列为 0，因此直接令自由分量取 `1..max_coefficient` 的整数，
主元分量按 `x_i = Σ b_k[i]·t_k` 用精确有理数计算，保留全部主元分量也是
`1..max_coefficient` 整数且整体 gcd 为 1 的向量。每个主元行预先通分（分母
LCM），热循环只剩整数乘加；自由变量最多 3 个、系数上限 30，枚举规模至多
`30³ = 27000` 组，通常在 0.2 秒内完成。

全程使用有理数，不会出现浮点误差，大系数反应（如
`10 K4Fe(CN)6 + 122 KMnO4 + 299 H2SO4 -> …`）也能精确求解。

---

## 5. 错误码

错误响应统一形如 `{"error": {"code": ..., "message": ..., "details": {...}}, "reaction_id": ..., "cached": ...}`。

| code | HTTP | 含义 |
| --- | --- | --- |
| `FORMULA_PARSE_ERROR` | 400 | 化学式无法解析（未知元素、括号不匹配、电荷非法等） |
| `DUPLICATE_SPECIES` | 409 | 同一侧出现重复物种 |
| `INVALID_REQUEST` | 422 | 请求体缺失/为空/物种过多（>100）；`/stoichiometry/batch` 投料缺失、重复、含未知物种或摩尔数非法（含 JSON 数字）；`/stoichiometry/plan` 目标重复、含非生成物物种或摩尔数非法 |
| `IMPOSSIBLE_BALANCE` | 422 | 守恒方程组无解（含电荷无法守恒） |
| `NON_UNIQUE_BALANCE` | 422 | 零空间维数 > 1，配平不唯一（返回 `nullity` 与自由度） |
| `ZERO_COEFFICIENT` | 422 | 唯一解要求某些物种系数为 0 |
| `INVALID_STOICHIOMETRY` | 422 | 唯一解符号不一致，物种可能需要换边 |
| `VARIANT_SPACE_TOO_LARGE` | 422 | `/balance/variants` 零空间维数（自由变量数）> 3，枚举空间过大（返回 `nullity`） |
| `NOT_FOUND` | 404 | `reaction_id` 不在缓存中 |

不唯一示例：`C + O2 -> CO + CO2` 缺少约束，返回

```json
{
  "error": {
    "code": "NON_UNIQUE_BALANCE",
    "message": "配平不唯一：零空间维数为 2，存在 1 个自由度……",
    "details": {"nullity": 2, "degrees_of_freedom": 1, "elements": ["C", "O"]}
  }
}
```

---

## 6. 项目结构

```
app/
  parser.py         递归下降化学式解析（元素/括号/电荷/电子）
  balancer.py       守恒矩阵、有理数 RREF、零空间、最小正整数、多变体枚举、守恒视图
  cache.py          规范化签名（SHA-256）与 LRU 幂等缓存
  stoichiometry.py  投料量 -> 反应进度（正向批量计量）与目标产量 -> 最小反应进度/理论投料（反向规划），Decimal 精确计算
  service.py        业务编排与响应组装
  models.py         Pydantic 请求/响应模型
  main.py           FastAPI 路由与异常处理
tests/              解析器、配平核心、HTTP API、批量计量、反向投料规划共 175 个测试
Dockerfile
docker-compose.yml
requirements.txt / requirements-dev.txt
```

## 7. 已验证的反应示例

```
2H2 + O2 -> 2H2O
4Fe + 3O2 -> 2Fe2O3
C3H8 + 5O2 -> 3CO2 + 4H2O
2KMnO4 + 16HCl -> 2KCl + 2MnCl2 + 5Cl2 + 8H2O
3Ca(OH)2 + 2H3PO4 -> Ca3(PO4)2 + 6H2O
(NH4)2Cr2O7 -> Cr2O3 + N2 + 4H2O
Fe^2+ -> Fe^3+ + e^-
MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O
3Cu + 2NO3^- + 8H^+ -> 3Cu^2+ + 2NO + 4H2O
10K4Fe(CN)6 + 122KMnO4 + 299H2SO4
  -> 162KHSO4 + 5Fe2(SO4)3 + 122MnSO4 + 60HNO3 + 60CO2 + 188H2O
```
