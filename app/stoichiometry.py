"""限量试剂与反应进度的精确计量计算。

给定唯一配平后的反应（最小互质正整数系数）与每种反应物的投料摩尔数
（:class:`decimal.Decimal`），计算：

* 反应进度 ξ = min(投料量 ÷ 配平系数)，系数为正整数，结果为非负小数；
* 所有并列取得最小值的限量试剂（可多于一种，如投料恰好按计量比）；
* 每种反应物的消耗量（ξ·系数）与剩余量（投料量 − 消耗量）；
* 每种生成物的理论产量（ξ·系数）。

为避免浮点误差，输入十进制先无损转成 :class:`fractions.Fraction`，
比较与最小值判定全部使用精确有理数；最终结果再按
:class:`decimal.Decimal`（28 位有效数字）格式化为字符串返回。
有限小数（分母仅含因子 2、5）给出精确值，无限小数按 28 位有效数字
银行家舍入。
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
from typing import Iterable

from .balancer import Species
from .cache import species_key

# 无限小数输出时保留的有效数字位数
SIGNIFICANT_DIGITS = 28

# 允许的数量级范围（Decimal 的 10 进制指数），防止 1e99999… 撑爆有理数内存
MAX_ABS_EXPONENT = 10_000


def parse_moles(text: str) -> Decimal:
    """把客户端的十进制字符串解析为非负 :class:`Decimal`。

    拒绝非法输入、NaN / Infinity、负数与过大数量级；允许千位下划线
    与科学计数法（``Decimal`` 原生语法）。
    """
    value = Decimal(text.strip())
    if not value.is_finite():
        raise ValueError("摩尔数不能是 NaN 或无穷大")
    # is_signed() 同时覆盖负数与 "-0" / "-0.0" 这类负零写法
    if value.is_signed() and value <= 0:
        raise ValueError("摩尔数不能为负数")
    if value.adjusted() > MAX_ABS_EXPONENT or value.as_tuple().exponent < -MAX_ABS_EXPONENT:
        raise ValueError(f"摩尔数的数量级超出允许范围（±1e{MAX_ABS_EXPONENT}）")
    return value


def decimal_to_fraction(value: Decimal) -> Fraction:
    """Decimal -> Fraction 的无损精确转换（不经过 float）。"""
    sign, digits, exponent = value.as_tuple()
    numerator = int("".join(map(str, digits)))
    if sign:
        numerator = -numerator
    if exponent >= 0:
        return Fraction(numerator * 10**exponent, 1)
    return Fraction(numerator, 10**-exponent)


def _factor_count(n: int, p: int) -> int:
    count = 0
    while n % p == 0:
        n //= p
        count += 1
    return count


def _terminating_places(denominator: int) -> int:
    """分母只含 2、5 时，精确十进制表示所需的小数位数 max(v2, v5)。"""
    twos = _factor_count(denominator, 2)
    rest = denominator // 2**twos
    fives = _factor_count(rest, 5)
    return max(twos, fives)


def _is_terminating(denominator: int) -> bool:
    twos = _factor_count(denominator, 2)
    rest = denominator // 2**twos
    rest //= 5 ** _factor_count(rest, 5)
    return rest == 1


def _fixed_point(decimal_value: Decimal) -> str:
    """把 Decimal 渲染为不带指数、不带尾随零的定点字符串。"""
    sign, digits, exponent = decimal_value.as_tuple()
    digit_text = "".join(map(str, digits))
    if exponent >= 0:
        body = digit_text + "0" * exponent
    else:
        places = -exponent
        if places < len(digit_text):
            body = digit_text[:-places] + "." + digit_text[-places:]
        else:
            body = "0." + "0" * (places - len(digit_text)) + digit_text
        if "." in body:
            body = body.rstrip("0").rstrip(".")
    return ("-" if sign and digit_text != "0" else "") + (body or "0")


def fraction_to_decimal_string(value: Fraction) -> str:
    """把精确有理数量渲染为十进制字符串。

    分母约分后只含 2、5 因子时输出精确有限小数；否则按
    :data:`SIGNIFICANT_DIGITS` 位有效数字舍入（ROUND_HALF_EVEN）。
    始终使用定点记法（不出现科学计数法），并去掉多余的尾随零。
    """
    if value == 0:
        return "0"
    denominator = value.denominator
    terminating = _is_terminating(denominator)
    places = _terminating_places(denominator) if terminating else 0

    with localcontext() as ctx:
        if terminating:
            # 精确除法：精度覆盖整数部分与所需小数位
            integer_digits = max(len(str(value.numerator)) - len(str(denominator)) + 1, 1)
            ctx.prec = integer_digits + places + 2
            decimal_value = Decimal(value.numerator) / Decimal(denominator)
            if places:
                decimal_value = decimal_value.quantize(Decimal(1).scaleb(-places))
        else:
            ctx.prec = SIGNIFICANT_DIGITS
            decimal_value = Decimal(value.numerator) / Decimal(denominator)
    return _fixed_point(decimal_value)


def compute_stoichiometry(
    reactants: Iterable[Species],
    products: Iterable[Species],
    coefficients: dict[str, int],
    feed: dict[str, Decimal],
) -> dict:
    """计算反应进度、限量试剂、各反应物消耗/剩余与各生成物理论产量。

    ``coefficients`` 与 ``feed`` 均以 :func:`app.cache.species_key` 为键；
    ``feed`` 必须覆盖全部反应物且值非负。

    物种按规范键排序输出，因此同一反应调换书写顺序得到的结果完全一致。
    """
    reactants = sorted(reactants, key=lambda sp: sp.formula.key)
    products = sorted(products, key=lambda sp: sp.formula.key)

    feed_fractions = {key: decimal_to_fraction(amount) for key, amount in feed.items()}

    # 反应进度 ξ = min(投料量 / 系数)
    ratios: dict[str, Fraction] = {
        species_key(sp): feed_fractions[species_key(sp)]
        / coefficients[species_key(sp)]
        for sp in reactants
    }
    extent = min(ratios.values())
    limiting = {key for key, ratio in ratios.items() if ratio == extent}

    reactant_blocks: list[dict] = []
    for sp in reactants:
        key = species_key(sp)
        coef = coefficients[key]
        fed = feed_fractions[key]
        consumed = extent * coef
        remaining = fed - consumed
        reactant_blocks.append(
            {
                "input": sp.raw,
                "normalized": sp.normalized,
                "coefficient": coef,
                "fed": fraction_to_decimal_string(fed),
                "consumed": fraction_to_decimal_string(consumed),
                "remaining": fraction_to_decimal_string(remaining),
                "is_limiting": key in limiting,
            }
        )

    product_blocks: list[dict] = []
    for sp in products:
        key = species_key(sp)
        coef = coefficients[key]
        produced = extent * coef
        product_blocks.append(
            {
                "input": sp.raw,
                "normalized": sp.normalized,
                "coefficient": coef,
                "theoretical_yield": fraction_to_decimal_string(produced),
            }
        )

    return {
        "reaction_extent": {
            "value": fraction_to_decimal_string(extent),
            "unit": "mol",
        },
        "limiting_reactants": [block["normalized"] for block in reactant_blocks
                               if block["is_limiting"]],
        "reactants": reactant_blocks,
        "products": product_blocks,
    }
