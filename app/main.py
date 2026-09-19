"""FastAPI 应用入口：环状 DNA 双酶切鉴别服务。

核心功能为纯本地计算，不依赖数据库、缓存或任何外部服务。
"""
from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from app.config import settings
from app.digest import recommend
from app.identify import identify
from app.models import (
    DigestRequest,
    DigestResponse,
    HealthResponse,
    IdentifyRequest,
    IdentifyResponse,
)
from app.validation import (
    DigestValidationError,
    validate_identify,
    validate_request,
)

VERSION = "1.0.0"

app = FastAPI(
    title="Circular DNA Double-Digest Discrimination Service",
    description=(
        "圆形（环状）DNA 双酶切鉴别服务：枚举不同酶对，按圆环相邻切点计算"
        "参考型与变异型的片段长度，按分辨容差聚类条带，并推荐可区分两种样本的酶对。"
    ),
    version=VERSION,
)


def _validation_response(errors: list[dict]) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={"code": "validation_failed", "errors": errors},
    )


@app.exception_handler(DigestValidationError)
async def digest_validation_handler(request: Request, exc: DigestValidationError) -> JSONResponse:
    return _validation_response(exc.detail["errors"])


@app.exception_handler(RequestValidationError)
async def request_validation_handler(request: Request, exc: RequestValidationError) -> JSONResponse:
    errors = [
        {"loc": [str(p) for p in err["loc"]], "msg": err["msg"], "type": err["type"]}
        for err in exc.errors()
    ]
    return _validation_response(errors)


@app.get("/health", response_model=HealthResponse, tags=["health"], summary="健康检查")
async def health() -> HealthResponse:
    return HealthResponse(status="ok", service=settings.APP_NAME, version=VERSION)


@app.post(
    "/api/v1/digest/recommend",
    response_model=DigestResponse,
    tags=["digest"],
    summary="枚举酶对并推荐可区分参考型与变异型的双酶切方案",
)
async def digest_recommend(req: DigestRequest) -> DigestResponse:
    validate_request(req)
    return recommend(req)


@app.post(
    "/api/v1/digest/identify",
    response_model=IdentifyResponse,
    tags=["digest"],
    summary="凝胶观测条带与参考型/变异型预测条带的一对一匹配鉴别",
)
async def digest_identify(req: IdentifyRequest) -> IdentifyResponse:
    validate_identify(req)
    return identify(req)


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("app.main:app", host=settings.HOST, port=settings.PORT)
