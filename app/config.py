"""运行期配置：全部通过环境变量注入，均有默认值，核心功能不依赖任何外部服务。"""
import os


class Settings:
    APP_NAME: str = os.getenv("APP_NAME", "circular-dna-digest")
    # 服务监听地址与端口（README / Compose 中使用 8000）
    HOST: str = os.getenv("APP_HOST", "0.0.0.0")
    PORT: int = int(os.getenv("APP_PORT", "8000"))

    # 单次请求允许提交的最大酶数量，避免组合枚举过大
    MAX_ENZYMES: int = int(os.getenv("MAX_ENZYMES", "64"))


settings = Settings()
