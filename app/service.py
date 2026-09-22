"""业务编排层：解析、配平、缓存与响应组装。"""

from __future__ import annotations

from decimal import Decimal, InvalidOperation
from typing import Any

from .balancer import (
    VARIANT_SPACE_TOO_LARGE,
    BalanceError,
    Species,
    build_balanced_view,
    build_raw_view,
    enumerate_variants_parsed,
    parse_side,
    solve_reaction,
)
from .cache import ReactionRecord, build_signature, cache, reaction_id_for, species_key
from .parser import FormulaError, parse_formula
from .stoichiometry import compute_reverse_plan, compute_stoichiometry, parse_moles


# 错误码 -> HTTP 状态码
STATUS_CODES = {
    "FORMULA_PARSE_ERROR": 400,
    "DUPLICATE_SPECIES": 409,
    "INVALID_REQUEST": 422,
    "IMPOSSIBLE_BALANCE": 422,
    "NON_UNIQUE_BALANCE": 422,
    "ZERO_COEFFICIENT": 422,
    "INVALID_STOICHIOMETRY": 422,
    "CHARGE_NOT_CONSERVED": 422,
    VARIANT_SPACE_TOO_LARGE: 422,
}


def error_payload(code: str, message: str, details: dict[str, Any]) -> dict[str, Any]:
    return {"code": code, "message": message, "details": details}


def _parse_all(
    reactant_texts: list[str], product_texts: list[str]
) -> tuple[list[Species], list[Species]]:
    try:
        return list(parse_side(reactant_texts, "反应物")), list(
            parse_side(product_texts, "生成物")
        )
    except FormulaError as exc:
        raise BalanceError("FORMULA_PARSE_ERROR", str(exc)) from exc


def _equation(
    reactants: list[Species],
    products: list[Species],
    coefficient_map: dict[str, int] | None,
) -> str:
    def term(sp: Species) -> str:
        coef = coefficient_map.get(species_key(sp), 1) if coefficient_map else 1
        return ("" if coef == 1 else str(coef)) + sp.normalized

    return " + ".join(map(term, reactants)) + " -> " + " + ".join(map(term, products))


def _normalized_block(
    reactants: list[Species],
    products: list[Species],
    coefficient_map: dict[str, int] | None,
) -> dict[str, Any]:
    return {
        "reactants": [sp.normalized for sp in reactants],
        "products": [sp.normalized for sp in products],
        "equation": _equation(reactants, products, coefficient_map),
    }


def _species_coefficients(
    side: list[Species], coefficient_map: dict[str, int]
) -> list[dict[str, Any]]:
    return [
        {
            "input": sp.raw,
            "normalized": sp.normalized,
            "coefficient": coefficient_map[species_key(sp)],
        }
        for sp in side
    ]


def _coef_map(solved) -> dict[str, int]:
    return {
        species_key(sp): coef
        for sp, coef in zip(solved.all_species, solved.coefficients)
    }


# ------------------------------------------------------------ 解析失败处理
def _text_signature(reactant_texts: list[str], product_texts: list[str]) -> str:
    """解析失败时的退化幂等键：压缩空白后排序。"""
    def norm(texts: list[str]) -> list[str]:
        return sorted("".join(t.split()) for t in texts)

    return f"raw:r={norm(reactant_texts)!r};p={norm(product_texts)!r}"


def _record_parse_error(
    reactant_texts: list[str], product_texts: list[str], exc: BalanceError
) -> tuple[ReactionRecord, bool]:
    signature = _text_signature(reactant_texts, product_texts)
    existing = cache.get_by_signature(signature)
    if existing is not None:
        return existing, True
    record = cache.put(
        signature, [], [], None, [],
        error_payload(exc.code, exc.message, dict(exc.details)),
    )
    return record, False


# ------------------------------------------------------------ 配平主流程
def process_balance(
    reactant_texts: list[str], product_texts: list[str]
) -> tuple[int, dict[str, Any]]:
    """配平入口，返回 (http_status, payload)。成功与错误结果都按规范签名缓存。"""
    try:
        reactants, products = _parse_all(reactant_texts, product_texts)
    except BalanceError as exc:
        record, cached = _record_parse_error(reactant_texts, product_texts, exc)
        return STATUS_CODES[record.error["code"]], {
            "error": record.error,
            "reaction_id": record.reaction_id,
            "cached": cached,
        }

    signature = build_signature(reactants, products)
    record = cache.get_by_signature(signature)
    cached = record is not None

    if record is None:
        try:
            solved = solve_reaction(reactant_texts, product_texts)
            record = cache.put(
                signature,
                list(solved.reactants),
                list(solved.products),
                _coef_map(solved),
                list(solved.elements),
                None,
            )
        except BalanceError as exc:
            record = cache.put(
                signature,
                reactants,
                products,
                None,
                list(exc.details.get("elements", [])),
                error_payload(exc.code, exc.message, dict(exc.details)),
            )

    if record.error is not None:
        return STATUS_CODES[record.error["code"]], {
            "error": record.error,
            "reaction_id": record.reaction_id,
            "cached": cached,
        }
    return 200, build_balance_payload(record, reactants, products, cached)


def build_balance_payload(
    record: ReactionRecord,
    reactants: list[Species],
    products: list[Species],
    cached: bool,
) -> dict[str, Any]:
    coefficient_map = record.coefficients
    charge_left = sum(
        coefficient_map[species_key(sp)] * sp.formula.charge for sp in reactants
    )
    charge_right = sum(
        coefficient_map[species_key(sp)] * sp.formula.charge for sp in products
    )
    return {
        "reaction_id": record.reaction_id,
        "cached": cached,
        "normalized": _normalized_block(reactants, products, coefficient_map),
        "coefficients": {
            "reactants": _species_coefficients(reactants, coefficient_map),
            "products": _species_coefficients(products, coefficient_map),
            "vector": [
                coefficient_map[species_key(sp)] for sp in (*reactants, *products)
            ],
        },
        "balanced": True,
        "elements": list(record.elements),
        "charge_balance": {
            "reactant_total": charge_left,
            "product_total": charge_right,
            "difference": charge_left - charge_right,
            "balanced": charge_left == charge_right,
        },
    }


# ------------------------------------------------------------ 逐元素分析
def process_analyze(
    reactant_texts: list[str], product_texts: list[str]
) -> tuple[int, dict[str, Any]]:
    """逐元素差异分析。可解析即返回 200（含原始/平衡两个视图）；无法解析返回 400。"""
    try:
        reactants, products = _parse_all(reactant_texts, product_texts)
    except BalanceError as exc:
        record, cached = _record_parse_error(reactant_texts, product_texts, exc)
        return STATUS_CODES[record.error["code"]], {
            "reaction_id": record.reaction_id,
            "cached": cached,
            "error": record.error,
        }

    signature = build_signature(reactants, products)
    record = cache.get_by_signature(signature)
    cached = record is not None
    elements = sorted(
        {el for sp in (*reactants, *products) for el in sp.formula.composition}
    )

    if record is None:
        try:
            solved = solve_reaction(reactant_texts, product_texts)
            record = cache.put(
                signature,
                list(solved.reactants),
                list(solved.products),
                _coef_map(solved),
                list(solved.elements),
                None,
            )
        except BalanceError as exc:
            record = cache.put(
                signature,
                reactants,
                products,
                None,
                list(exc.details.get("elements", [])),
                error_payload(exc.code, exc.message, dict(exc.details)),
            )

    raw = build_raw_view(reactants, products, elements)
    payload: dict[str, Any] = {
        "reaction_id": record.reaction_id,
        "cached": cached,
        "normalized": _normalized_block(
            reactants, products, None if record.error else record.coefficients
        ),
        "raw": _view(reactants, products, raw, coefficient_map=None, assumed=1),
        "balanced_ok": record.error is None,
        "error": record.error,
    }
    if record.error is None:
        solved = solve_reaction(reactant_texts, product_texts)
        balanced_view = build_balanced_view(solved)
        payload["balanced"] = _view(
            list(solved.reactants),
            list(solved.products),
            balanced_view,
            coefficient_map=record.coefficients,
            assumed=None,
        )
    else:
        payload["balanced"] = None
    return 200, payload


def _view(
    reactants: list[Species],
    products: list[Species],
    view: dict[str, Any],
    *,
    coefficient_map: dict[str, int] | None,
    assumed: int | None,
) -> dict[str, Any]:
    def side_block(side: list[Species]) -> list[dict[str, Any]]:
        block = []
        for sp in side:
            coef = (
                coefficient_map[species_key(sp)] if coefficient_map else (assumed or 1)
            )
            block.append(
                {"input": sp.raw, "normalized": sp.normalized, "coefficient": coef}
            )
        return block

    result: dict[str, Any] = {
        "species": {
            "reactants": side_block(reactants),
            "products": side_block(products),
        },
        "elements": view["elements"],
        "charge": view["charge"],
    }
    if assumed is not None:
        result["coefficients_assumed"] = assumed
    return result


def get_reaction(reaction_id: str) -> dict[str, Any] | None:
    """按 reaction_id 取回缓存结果（幂等查询）。"""
    record = cache.get(reaction_id)
    if record is None:
        return None
    reactants = _record_species(record, "canonical_reactants")
    products = _record_species(record, "canonical_products")
    if record.error is not None:
        return {
            "reaction_id": record.reaction_id,
            "cached": True,
            "error": record.error,
        }
    return build_balance_payload(record, reactants, products, True)


def _record_species(record: ReactionRecord, attr: str) -> list[Species]:
    """缓存里只存了 normalized 文本；查询回显时构造轻量 Species。"""
    from .parser import parse_formula

    return [
        Species(
            raw=item["normalized"],
            normalized=item["normalized"],
            formula=parse_formula(item["normalized"]),
        )
        for item in getattr(record, attr)
    ]


# ------------------------------------------------------------ 多变体枚举
def process_variants(
    reactant_texts: list[str],
    product_texts: list[str],
    max_coefficient: int,
    limit: int,
) -> tuple[int, dict[str, Any]]:
    """枚举配平变体，返回 (http_status, payload)。结果不写幂等缓存。"""
    try:
        reactants, products = _parse_all(reactant_texts, product_texts)
    except BalanceError as exc:
        # 解析 / 去重阶段失败：沿用与 /balance 一致的退化签名生成 ID
        reaction_id = reaction_id_for(
            _text_signature(reactant_texts, product_texts)
        )
        return STATUS_CODES[exc.code], {
            "error": error_payload(exc.code, exc.message, dict(exc.details)),
            "reaction_id": reaction_id,
        }

    try:
        variant_set = enumerate_variants_parsed(reactants, products, max_coefficient)
    except BalanceError as exc:
        reaction_id = reaction_id_for(build_signature(reactants, products))
        return STATUS_CODES[exc.code], {
            "error": error_payload(exc.code, exc.message, dict(exc.details)),
            "reaction_id": reaction_id,
        }

    signature = build_signature(reactants, products)
    reaction_id = reaction_id_for(signature)

    total = len(variant_set.coefficients)
    page = variant_set.coefficients[:limit]

    items: list[dict[str, Any]] = []
    for vec in page:
        coefficient_map = {
            species_key(sp): coef
            for sp, coef in zip((*reactants, *products), vec)
        }
        items.append(
            {
                "coefficients": {
                    "reactants": _species_coefficients(reactants, coefficient_map),
                    "products": _species_coefficients(products, coefficient_map),
                    "vector": list(vec),
                },
                "sum": sum(vec),
                "equation": _equation(reactants, products, coefficient_map),
            }
        )

    return 200, {
        "reaction_id": reaction_id,
        "normalized": _normalized_block(reactants, products, None),
        "max_coefficient": max_coefficient,
        "limit": limit,
        "nullity": variant_set.nullity,
        "count": len(page),
        "total": total,
        "truncated": total > limit,
        "variants": items,
    }


# ------------------------------------------------------------ 限量试剂批量计量
def process_stoichiometry(
    reactant_texts: list[str],
    product_texts: list[str],
    amount_items: list[dict[str, str]],
) -> tuple[int, dict[str, Any]]:
    """限量试剂批量计量，返回 (http_status, payload)。

    先复用 :func:`process_balance` 的唯一配平、规范化与幂等缓存逻辑；
    配平失败（解析错误、重复物种、不可唯一配平等）沿用其结构化错误响应。
    随后校验投料与规范化反应物一一对应，并以精确十进制计算反应进度、
    并列限量试剂、各反应物消耗/剩余与各生成物理论产量。
    """
    status, balance_payload = process_balance(reactant_texts, product_texts)
    if status != 200:
        return status, balance_payload

    reaction_id = balance_payload["reaction_id"]
    cached = balance_payload["cached"]

    # 规范化反应物（parse_side 已按规范键排序，与配平向量顺序一致）
    reactants = list(parse_side(reactant_texts, "反应物"))
    products = list(parse_side(product_texts, "生成物"))
    coefficient_map = {
        species_key(sp): coef
        for sp, coef in zip(
            (*reactants, *products), balance_payload["coefficients"]["vector"]
        )
    }
    reactant_keys = {species_key(sp) for sp in reactants}

    # ---- 投料条目：逐条解析物种与摩尔数，聚合所有可定位的问题
    parsed_entries: list[tuple[str, Species]] = []
    invalid_moles: list[dict[str, str]] = []
    feed: dict[str, Decimal] = {}
    duplicate: list[str] = []
    seen_keys: set[str] = set()

    for item in amount_items:
        species_text = item["species"]
        moles_text = item["moles"]
        try:
            formula = parse_formula(species_text)
        except FormulaError as exc:
            # 投料物种本身无法解析：沿用化学式解析错误（400）
            return STATUS_CODES["FORMULA_PARSE_ERROR"], {
                "reaction_id": reaction_id,
                "cached": cached,
                "error": error_payload(
                    "FORMULA_PARSE_ERROR",
                    str(exc),
                    {"field": "amounts.species", "input": species_text},
                ),
            }
        sp = Species(raw=species_text, normalized=formula.normalized, formula=formula)
        key = species_key(sp)

        if key in seen_keys:
            duplicate.append(formula.normalized)
        seen_keys.add(key)

        try:
            amount = parse_moles(moles_text)
        except (ValueError, InvalidOperation):
            invalid_moles.append({"species": formula.normalized, "moles": moles_text})
        else:
            parsed_entries.append((moles_text, sp))
            feed.setdefault(key, amount)

    unknown = sorted(
        {sp.normalized for _, sp in parsed_entries if species_key(sp) not in reactant_keys}
    )
    provided_keys = {species_key(sp) for _, sp in parsed_entries}
    normalized_by_key = {species_key(sp): sp.normalized for sp in reactants}
    missing = [
        normalized_by_key[key] for key in reactant_keys if key not in provided_keys
    ]

    if duplicate or unknown or missing or invalid_moles:
        details: dict[str, Any] = {}
        if missing:
            details["missing"] = sorted(missing)
        if duplicate:
            details["duplicate"] = sorted(set(duplicate))
        if unknown:
            details["unknown"] = unknown
        if invalid_moles:
            details["invalid_moles"] = invalid_moles
        return STATUS_CODES["INVALID_REQUEST"], {
            "reaction_id": reaction_id,
            "cached": cached,
            "error": error_payload(
                "INVALID_REQUEST",
                "投料记录必须与规范化反应物一一对应：不能有缺失、重复或未知物种，"
                "且每条摩尔数必须是非负十进制字符串（不能为 JSON 数字）",
                details,
            ),
        }

    result = compute_stoichiometry(reactants, products, coefficient_map, feed)
    return 200, {
        "reaction_id": reaction_id,
        "cached": cached,
        "normalized": balance_payload["normalized"],
        **result,
    }


# ------------------------------------------------------------ 反向投料规划
def process_reverse_plan(
    reactant_texts: list[str],
    product_texts: list[str],
    target_items: list[dict[str, str]],
) -> tuple[int, dict[str, Any]]:
    """反向投料规划，返回 (http_status, payload)。

    先复用 :func:`process_balance` 的唯一配平、规范化与幂等缓存逻辑；
    配平失败沿用其结构化错误响应。随后校验目标条目都能映射到规范化
    生成物（允许只提交部分生成物，但不能重复或出现未知物种），并以精确
    十进制计算最小反应进度 ξ = max(目标量 ÷ 系数)、各反应物理论需求量、
    各目标产物达成量/超额量与全部生成物达成量。
    """
    status, balance_payload = process_balance(reactant_texts, product_texts)
    if status != 200:
        return status, balance_payload

    reaction_id = balance_payload["reaction_id"]
    cached = balance_payload["cached"]

    # parse_side 已按规范键排序，与配平向量顺序一致
    reactants = list(parse_side(reactant_texts, "反应物"))
    products = list(parse_side(product_texts, "生成物"))
    coefficient_map = {
        species_key(sp): coef
        for sp, coef in zip(
            (*reactants, *products), balance_payload["coefficients"]["vector"]
        )
    }
    product_keys = {species_key(sp) for sp in products}

    # ---- 目标条目：逐条解析物种与摩尔数，聚合所有可定位的问题
    parsed_entries: list[tuple[Species, Decimal]] = []
    invalid_moles: list[dict[str, str]] = []
    duplicate: list[str] = []
    seen_keys: set[str] = set()

    for item in target_items:
        species_text = item["species"]
        moles_text = item["moles"]
        try:
            formula = parse_formula(species_text)
        except FormulaError as exc:
            # 目标物种本身无法解析：沿用化学式解析错误（400）
            return STATUS_CODES["FORMULA_PARSE_ERROR"], {
                "reaction_id": reaction_id,
                "cached": cached,
                "error": error_payload(
                    "FORMULA_PARSE_ERROR",
                    str(exc),
                    {"field": "targets.species", "input": species_text},
                ),
            }
        sp = Species(raw=species_text, normalized=formula.normalized, formula=formula)
        key = species_key(sp)

        if key in seen_keys:
            duplicate.append(formula.normalized)
        seen_keys.add(key)

        try:
            amount = parse_moles(moles_text)
        except (ValueError, InvalidOperation):
            invalid_moles.append({"species": formula.normalized, "moles": moles_text})
        else:
            parsed_entries.append((sp, amount))

    unknown = sorted(
        {sp.normalized for sp, _ in parsed_entries if species_key(sp) not in product_keys}
    )

    if duplicate or unknown or invalid_moles:
        details: dict[str, Any] = {}
        if duplicate:
            details["duplicate"] = sorted(set(duplicate))
        if unknown:
            details["unknown"] = unknown
        if invalid_moles:
            details["invalid_moles"] = invalid_moles
        return STATUS_CODES["INVALID_REQUEST"], {
            "reaction_id": reaction_id,
            "cached": cached,
            "error": error_payload(
                "INVALID_REQUEST",
                "目标产物必须是规范化生成物之一：不能重复或包含未知物种"
                "（允许只提交部分生成物），且每条目标摩尔数必须是非负十进制"
                "字符串（不能为 JSON 数字）",
                details,
            ),
        }

    result = compute_reverse_plan(
        reactants, products, coefficient_map, parsed_entries
    )
    return 200, {
        "reaction_id": reaction_id,
        "cached": cached,
        "normalized": balance_payload["normalized"],
        **result,
    }
