"""FastAPI 应用入口。

路由（统一前缀 /api/v1）：

* ``POST /balance``          配平并返回最小正整数系数；
* ``POST /balance/variants`` 枚举零空间内全部互质正整数配平方案；
* ``POST /analyze``          逐元素差异说明（配平前 + 配平后视图）；
* ``POST /stoichiometry/batch`` 投料摩尔数 -> 反应进度、并列限量试剂、
  各反应物消耗/剩余与各生成物理论产量；
* ``POST /stoichiometry/plan``  生成物目标摩尔数 -> 最小反应进度、
  各反应物理论需求量、各目标产物达成量与超额量；
* ``GET  /reactions/{reaction_id}``  按幂等 ID 取回结果；
* ``GET  /health``   健康检查。
"""

from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from . import __version__
from .models import (
    AnalyzeResponse,
    BalanceResponse,
    ReactionRequest,
    ReversePlanRequest,
    ReversePlanResponse,
    StoichiometryResponse,
    StoichiometryRequest,
    VariantRequest,
    VariantsResponse,
)
from .service import (
    STATUS_CODES,
    get_reaction,
    process_analyze,
    process_balance,
    process_reverse_plan,
    process_stoichiometry,
    process_variants,
)

app = FastAPI(
    title="化学反应计量平衡服务",
    description=(
        "解析含元素、括号与电荷的化学式，建立元素守恒与电荷守恒矩阵，"
        "通过有理数高斯-若尔当消元求解最小互质正整数系数。"
    ),
    version=__version__,
)


@app.exception_handler(RequestValidationError)
async def handle_validation_error(request: Request, exc: RequestValidationError):
    return JSONResponse(
        status_code=422,
        content={
            "error": {
                "code": "INVALID_REQUEST",
                "message": "请求体校验失败：请提供非空的 reactants 与 products 字符串数组",
                "details": {"validation": exc.errors()},
            }
        },
    )


@app.get("/health", tags=["系统"])
async def health():
    return {"status": "ok", "version": __version__}


@app.post(
    "/api/v1/balance",
    response_model=BalanceResponse,
    tags=["配平"],
    summary="配平化学反应",
    responses={400: {"description": "化学式解析失败"}, 409: {"description": "同侧重复物种"},
               422: {"description": "无法配平 / 请求非法"}},
)
async def balance(request: ReactionRequest):
    status, payload = process_balance(request.reactants, request.products)
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/balance/variants",
    response_model=VariantsResponse,
    tags=["配平"],
    summary="枚举所有互质正整数配平方案（自由变量不超过 3 个）",
    responses={400: {"description": "化学式解析失败"}, 409: {"description": "同侧重复物种"},
               422: {"description": "无法配平 / 自由变量超过 3 个 / 请求非法"}},
)
async def balance_variants(request: VariantRequest):
    status, payload = process_variants(
        request.reactants,
        request.products,
        request.max_coefficient,
        request.limit,
    )
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/analyze",
    response_model=AnalyzeResponse,
    tags=["配平"],
    summary="逐元素守恒差异分析",
    responses={400: {"description": "化学式解析失败"}},
)
async def analyze(request: ReactionRequest):
    status, payload = process_analyze(request.reactants, request.products)
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/stoichiometry/batch",
    response_model=StoichiometryResponse,
    tags=["计量"],
    summary="批量限量试剂计量：反应进度、消耗/剩余与理论产量",
    responses={
        400: {"description": "化学式解析失败"},
        409: {"description": "同侧重复物种 / 投料重复"},
        422: {"description": "无法唯一配平 / 投料与反应物不对应 / 摩尔数非法"},
    },
)
async def stoichiometry_batch(request: StoichiometryRequest):
    status, payload = process_stoichiometry(
        request.reactants,
        request.products,
        [item.model_dump() for item in request.amounts],
    )
    return JSONResponse(status_code=status, content=payload)


@app.post(
    "/api/v1/stoichiometry/plan",
    response_model=ReversePlanResponse,
    tags=["计量"],
    summary="反向投料规划：按生成物目标摩尔数反推最小反应进度与理论投料需求",
    responses={
        400: {"description": "化学式解析失败"},
        409: {"description": "同侧重复物种"},
        422: {"description": "无法唯一配平 / 目标与生成物不对应 / 摩尔数非法"},
    },
)
async def stoichiometry_plan(request: ReversePlanRequest):
    status, payload = process_reverse_plan(
        request.reactants,
        request.products,
        [item.model_dump() for item in request.targets],
    )
    return JSONResponse(status_code=status, content=payload)


@app.get(
    "/api/v1/reactions/{reaction_id}",
    tags=["配平"],
    summary="按反应 ID 查询缓存结果（幂等）",
    responses={404: {"description": "未找到该反应"}},
)
async def get_one(reaction_id: str):
    payload = get_reaction(reaction_id)
    if payload is None:
        return JSONResponse(
            status_code=404,
            content={
                "error": {
                    "code": "NOT_FOUND",
                    "message": f"未找到反应 {reaction_id}（缓存可能已淘汰）",
                    "details": {"reaction_id": reaction_id},
                }
            },
        )
    status = 200 if "error" not in payload else STATUS_CODES[payload["error"]["code"]]
    return JSONResponse(status_code=status, content=payload)
