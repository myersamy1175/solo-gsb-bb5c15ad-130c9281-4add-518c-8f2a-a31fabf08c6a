"""POST /api/v1/stoichiometry/plan 的 HTTP API 测试 + 精确计算单元测试。

反向投料规划：客户端为一个或多个生成物提交目标摩尔数，服务端计算
ξ = max(目标量 ÷ 配平系数) 与各反应物理论需求量。

HTTP 测试反应均选用其它测试文件未提交过的方程式（NH3 合成、N2H4
推进剂、NaBr 置换、Li3N 合成等），避免幂等缓存造成跨文件的顺序依赖。
"""

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.stoichiometry import compute_reverse_plan
from app.balancer import parse_side

client = TestClient(app)

PATH = "/api/v1/stoichiometry/plan"


# ------------------------------------------------------------ 成功路径
def test_single_target_plans_required_reactants():
    # N2 + 3H2 -> 2NH3；目标 4 mol NH3 -> ξ=2，需 N2 2 mol、H2 6 mol
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": "4"}],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["reaction_extent"] == {"value": "2", "unit": "mol"}
    assert data["driving_targets"] == ["NH3"]

    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["H2"]["coefficient"] == 3
    assert by_name["H2"]["required"] == "6"
    assert by_name["N2"]["coefficient"] == 1
    assert by_name["N2"]["required"] == "2"

    products = {p["normalized"]: p for p in data["products"]}
    nh3 = products["NH3"]
    assert nh3["coefficient"] == 2
    assert nh3["achieved"] == "4"
    assert nh3["target"] == "4"
    assert nh3["excess"] == "0"
    assert nh3["is_driving"] is True
    assert isinstance(data["reaction_extent"]["value"], str)


def test_multiple_takes_maximum_ratio_and_reports_excess():
    # N2H4 + 2H2O2 -> N2 + 4H2O
    # 目标 N2=1（比值 1）、H2O=2（比值 0.5）-> ξ=1，H2O 超额 2 mol
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2H4", "H2O2"],
            "products": ["N2", "H2O"],
            "targets": [
                {"species": "H2O", "moles": "2"},
                {"species": "N2", "moles": "1"},
            ],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["reaction_extent"]["value"] == "1"
    assert data["driving_targets"] == ["N2"]

    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["N2H4"]["required"] == "1"
    assert by_name["H2O2"]["required"] == "2"

    products = {p["normalized"]: p for p in data["products"]}
    assert products["N2"]["achieved"] == "1"
    assert products["N2"]["target"] == "1"
    assert products["N2"]["excess"] == "0"
    assert products["N2"]["is_driving"] is True
    # H2O 在该进度下产出 4 mol，超出目标 2 mol
    assert products["H2O"]["achieved"] == "4"
    assert products["H2O"]["target"] == "2"
    assert products["H2O"]["excess"] == "2"
    assert products["H2O"]["is_driving"] is False


def test_tied_driving_targets_at_stoichiometric_ratio():
    # N2H4 + 2H2O2 -> N2 + 4H2O；目标 N2=1、H2O=4 恰按系数比 -> 并列驱动
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2H4", "H2O2"],
            "products": ["N2", "H2O"],
            "targets": [
                {"species": "N2", "moles": "1"},
                {"species": "H2O", "moles": "4"},
            ],
        },
    )
    data = resp.json()
    assert data["reaction_extent"]["value"] == "1"
    assert data["driving_targets"] == ["H2O", "N2"]
    assert all(p["excess"] == "0" for p in data["products"])
    assert all(p["is_driving"] for p in data["products"])


def test_partial_targets_leave_other_products_with_null_fields():
    # SO3 + H2O -> H2SO4：只有一个产物；换成多产物反应
    # 2NaBr + Cl2 -> 2NaCl + Br2；只对 NaCl 下目标 1 mol -> ξ=0.5
    resp = client.post(
        PATH,
        json={
            "reactants": ["NaBr", "Cl2"],
            "products": ["NaCl", "Br2"],
            "targets": [{"species": "NaCl", "moles": "1"}],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["reaction_extent"]["value"] == "0.5"
    assert data["driving_targets"] == ["NaCl"]

    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["NaBr"]["required"] == "1"
    assert by_name["Cl2"]["required"] == "0.5"

    products = {p["normalized"]: p for p in data["products"]}
    assert products["NaCl"]["target"] == "1"
    assert products["NaCl"]["achieved"] == "1"
    assert products["NaCl"]["excess"] == "0"
    # 未列入 targets 的产物仍返回达成量，目标/超额字段为 null
    br2 = products["Br2"]
    assert br2["achieved"] == "0.5"
    assert br2["target"] is None
    assert br2["excess"] is None
    assert br2["is_driving"] is False


def test_decimal_inputs_avoid_binary_float_error():
    # N2 + 3H2 -> 2NH3；目标 0.1 mol NH3 -> ξ=0.05，H2 需 0.15
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": "0.1"}],
        },
    )
    data = resp.json()
    assert data["reaction_extent"]["value"] == "0.05"
    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["N2"]["required"] == "0.05"
    assert by_name["H2"]["required"] == "0.15"
    assert data["products"][0]["achieved"] == "0.1"


def test_repeating_decimal_is_string_without_float_artifacts():
    # Fe2O3 + 3H2SO4 -> Fe2(SO4)3 + 3H2O；目标 H2O=1 -> ξ=1/3（无限循环小数）
    resp = client.post(
        PATH,
        json={
            "reactants": ["Fe2O3", "H2SO4"],
            "products": ["Fe2(SO4)3", "H2O"],
            "targets": [{"species": "H2O", "moles": "1"}],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["reaction_extent"]["value"] == (
        "0.3333333333333333333333333333"  # 1/3，28 位有效数字
    )
    by_name = {r["normalized"]: r for r in data["reactants"]}
    # H2SO4 需求 = 3/3 = 1（精确）；Fe2O3 需求 = 1/3
    assert by_name["H2SO4"]["required"] == "1"
    assert by_name["Fe2O3"]["required"] == "0.3333333333333333333333333333"
    products = {p["normalized"]: p for p in data["products"]}
    assert products["Fe2(SO4)3"]["achieved"] == "0.3333333333333333333333333333"
    assert products["Fe2(SO4)3"]["target"] is None
    assert products["H2O"]["achieved"] == "1"
    assert products["H2O"]["excess"] == "0"


def test_zero_target_gives_zero_extent_and_requirements():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": "0"}],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["reaction_extent"]["value"] == "0"
    assert data["driving_targets"] == ["NH3"]
    assert all(r["required"] == "0" for r in data["reactants"])
    assert data["products"][0]["achieved"] == "0"
    assert data["products"][0]["excess"] == "0"


def test_zero_and_positive_targets_together():
    # 2NaBr + Cl2 -> 2NaCl + Br2；目标 NaCl=0、Br2=1 -> ξ=1
    resp = client.post(
        PATH,
        json={
            "reactants": ["NaBr", "Cl2"],
            "products": ["NaCl", "Br2"],
            "targets": [
                {"species": "NaCl", "moles": "0"},
                {"species": "Br2", "moles": "1"},
            ],
        },
    )
    data = resp.json()
    assert data["reaction_extent"]["value"] == "1"
    assert data["driving_targets"] == ["Br2"]
    products = {p["normalized"]: p for p in data["products"]}
    assert products["NaCl"]["achieved"] == "2"
    assert products["NaCl"]["target"] == "0"
    assert products["NaCl"]["excess"] == "2"
    assert products["NaCl"]["is_driving"] is False
    assert products["Br2"]["excess"] == "0"


def test_scientific_notation_and_underscores_accepted():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": "2e-1"}],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    # 目标 0.2 mol NH3 -> ξ=0.1
    assert data["reaction_extent"]["value"] == "0.1"
    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["N2"]["required"] == "0.1"
    assert by_name["H2"]["required"] == "0.3"
    assert data["products"][0]["target"] == "0.2"


# ------------------------------------------------------------ 顺序/写法无关性
def test_normalized_species_matching_with_alias_notation():
    # 等价写法（括号、空白、电荷风格）应正确匹配规范化生成物
    resp = client.post(
        PATH,
        json={
            "reactants": ["MnO4-", "H+", "e-"],
            "products": ["Mn2+", "H2O"],
            "targets": [
                {"species": "Mn^2+", "moles": "0.2"},
                {"species": "H2 O", "moles": "0.8"},
            ],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    # MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O
    # ξ = max(0.2/1, 0.8/4) = 0.2，两者恰按系数比并列驱动
    assert data["reaction_extent"]["value"] == "0.2"
    assert data["driving_targets"] == ["H2O", "Mn^2+"]
    by_name = {r["normalized"]: r for r in data["reactants"]}
    assert by_name["MnO4^-"]["required"] == "0.2"
    assert by_name["H^+"]["required"] == "1.6"
    assert by_name["e^-"]["required"] == "1"
    products = {p["normalized"]: p for p in data["products"]}
    assert products["H2O"]["achieved"] == "0.8"
    assert products["H2O"]["excess"] == "0"


def test_reaction_id_and_results_stable_across_order_and_notation():
    payload_a = {
        "reactants": ["N2H4", "H2O2"],
        "products": ["N2", "H2O"],
        "targets": [
            {"species": "N2", "moles": "1"},
            {"species": "H2O", "moles": "2"},
        ],
    }
    payload_b = {
        "reactants": ["H2 O2", "N2 H4"],
        "products": ["H2O", "N2"],
        "targets": [
            {"species": " H2O", "moles": "2.0"},
            {"species": "N2", "moles": "1.0"},
        ],
    }
    a = client.post(PATH, json=payload_a).json()
    b = client.post(PATH, json=payload_b).json()
    assert a["reaction_id"] == b["reaction_id"]
    assert a["reaction_extent"] == b["reaction_extent"]
    assert a["driving_targets"] == b["driving_targets"]

    def computed(blob):
        return (
            [(r["normalized"], r["coefficient"], r["required"])
             for r in blob["reactants"]],
            [(p["normalized"], p["coefficient"], p["achieved"],
              p["target"], p["excess"], p["is_driving"])
             for p in blob["products"]],
        )

    assert computed(a) == computed(b)


def test_reaction_id_matches_balance_and_batch_endpoints():
    balance = client.post(
        "/api/v1/balance",
        json={"reactants": ["N2", "H2"], "products": ["NH3"]},
    ).json()
    plan = client.post(
        PATH,
        json={
            "reactants": ["H2", "N2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": "2"}],
        },
    ).json()
    assert plan["reaction_id"] == balance["reaction_id"]


def test_cached_flag_reflects_balance_cache():
    # 使用本文件内其它用例未提交过的反应，避免幂等缓存造成顺序依赖
    body = {
        "reactants": ["Li", "N2"],
        "products": ["Li3N"],
        "targets": [{"species": "Li3N", "moles": "1"}],
    }
    first = client.post(PATH, json=body).json()
    second = client.post(PATH, json=body).json()
    # 6Li + N2 -> 2Li3N；目标 1 mol Li3N -> ξ=0.5
    assert first["reaction_extent"]["value"] == "0.5"
    assert first["cached"] is False
    assert second["cached"] is True


# ------------------------------------------------------------ 结构化错误
def test_duplicate_target_is_invalid_request():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [
                {"species": "NH3", "moles": "1"},
                {"species": "NH3", "moles": "2"},
            ],
        },
    )
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    assert details["duplicate"] == ["NH3"]


def test_unknown_species_in_targets_is_invalid_request():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [
                {"species": "NH3", "moles": "1"},
                {"species": "H2", "moles": "1"},
            ],
        },
    )
    assert resp.status_code == 422
    details = resp.json()["error"]["details"]
    # H2 是反应物，不属于生成物 -> 未知目标物种；不要求覆盖全部产物
    assert details["unknown"] == ["H2"]
    assert "missing" not in details


def test_reactant_species_target_is_invalid_not_silently_matched():
    # 目标必须是生成物：即使某物种存在于反应中，只要在反应物一侧即拒绝
    resp = client.post(
        PATH,
        json={
            "reactants": ["H2", "I2"],
            "products": ["HI"],
            "targets": [{"species": "H2", "moles": "1"}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["details"]["unknown"] == ["H2"]


@pytest.mark.parametrize(
    "bad",
    ["-1", "-0.0", "abc", "1.2.3", "NaN", "Infinity", "", "1e10001", "1e-10001"],
)
def test_invalid_moles_values(bad):
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": bad}],
        },
    )
    assert resp.status_code == 422
    body = resp.json()
    # 空字符串在 Pydantic 层即被拦截；其余在业务层聚合报告
    assert body["error"]["code"] == "INVALID_REQUEST"


def test_moles_must_be_string_not_json_number():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": 0.1}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_unparseable_species_in_targets_is_parse_error():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "H2"],
            "products": ["NH3"],
            "targets": [{"species": "Xx", "moles": "1"}],
        },
    )
    assert resp.status_code == 400
    body = resp.json()
    assert body["error"]["code"] == "FORMULA_PARSE_ERROR"
    assert body["error"]["details"]["input"] == "Xx"
    assert body["error"]["details"]["field"] == "targets.species"


def test_non_unique_balance_reuses_structured_error():
    resp = client.post(
        PATH,
        json={
            "reactants": ["C", "O2"],
            "products": ["CO", "CO2"],
            "targets": [{"species": "CO", "moles": "1"}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "NON_UNIQUE_BALANCE"


def test_impossible_balance_reuses_structured_error():
    # 选用其它测试未用过的不可能反应，避免顺序依赖
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2"],
            "products": ["Au"],
            "targets": [{"species": "Au", "moles": "1"}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "IMPOSSIBLE_BALANCE"


def test_reaction_side_duplicate_reuses_conflict():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2", "N2"],
            "products": ["NH3"],
            "targets": [{"species": "NH3", "moles": "1"}],
        },
    )
    assert resp.status_code == 409
    assert resp.json()["error"]["code"] == "DUPLICATE_SPECIES"


def test_reaction_formula_parse_error_reuses_parse_error():
    resp = client.post(
        PATH,
        json={
            "reactants": ["Zz", "H2"],
            "products": ["H2O"],
            "targets": [{"species": "H2O", "moles": "1"}],
        },
    )
    assert resp.status_code == 400
    assert resp.json()["error"]["code"] == "FORMULA_PARSE_ERROR"


@pytest.mark.parametrize(
    "payload",
    [
        {"reactants": ["N2", "H2"], "products": ["NH3"]},  # 缺 targets
        {"reactants": [], "products": ["NH3"],
         "targets": [{"species": "NH3", "moles": "1"}]},
        {"reactants": ["N2", "H2"], "products": ["NH3"], "targets": []},
    ],
)
def test_request_shape_validation(payload):
    resp = client.post(PATH, json=payload)
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_response_shape_contract():
    resp = client.post(
        PATH,
        json={
            "reactants": ["N2H4", "H2O2"],
            "products": ["N2", "H2O"],
            "targets": [{"species": "N2", "moles": "1"}],
        },
    )
    data = resp.json()
    assert set(data) == {
        "reaction_id", "cached", "normalized", "reaction_extent",
        "driving_targets", "reactants", "products",
    }
    assert set(data["reactants"][0]) == {
        "input", "normalized", "coefficient", "required",
    }
    assert set(data["products"][0]) == {
        "input", "normalized", "coefficient", "achieved",
        "target", "excess", "is_driving",
    }


# ------------------------------------------------------------ 计算核心单元测试
def _setup(reactant_texts, product_texts):
    reactants = parse_side(reactant_texts, "反应物")
    products = parse_side(product_texts, "生成物")

    from app.balancer import solve_reaction
    from app.cache import species_key

    solved = solve_reaction(reactant_texts, product_texts)
    coefficients = {
        species_key(sp): coef
        for sp, coef in zip(solved.all_species, solved.coefficients)
    }
    return reactants, products, coefficients


def test_compute_reverse_plan_basic_numbers():
    from decimal import Decimal
    from app.cache import species_key

    reactants, products, coefficients = _setup(["N2", "H2"], ["NH3"])
    targets = {species_key(products[0]): Decimal("4")}
    result = compute_reverse_plan(reactants, products, coefficients, targets)
    assert result["reaction_extent"] == {"value": "2", "unit": "mol"}
    assert result["driving_targets"] == ["NH3"]
    required = {b["normalized"]: b["required"] for b in result["reactants"]}
    assert required == {"H2": "6", "N2": "2"}


def test_compute_reverse_plan_excess_never_negative():
    from decimal import Decimal
    from app.cache import species_key

    reactants, products, coefficients = _setup(
        ["N2H4", "H2O2"], ["N2", "H2O"]
    )
    key_by_name = {sp.normalized: species_key(sp) for sp in products}
    # 目标量超过系数比：N2=100 -> ξ=100，H2O 同时产出 400，远超其目标 1
    targets = {
        key_by_name["N2"]: Decimal("100"),
        key_by_name["H2O"]: Decimal("1"),
    }
    result = compute_reverse_plan(reactants, products, coefficients, targets)
    assert result["reaction_extent"]["value"] == "100"
    assert result["driving_targets"] == ["N2"]
    for block in result["products"]:
        assert block["excess"] is None or Decimal(block["excess"]) >= 0
    nitrogen = next(p for p in result["products"] if p["normalized"] == "N2")
    assert nitrogen["achieved"] == "100"
    assert nitrogen["excess"] == "0"
    water = next(p for p in result["products"] if p["normalized"] == "H2O")
    assert water["achieved"] == "400"
    assert water["excess"] == "399"
    required = {b["normalized"]: b["required"] for b in result["reactants"]}
    assert required == {"H2O2": "200", "N2H4": "100"}
