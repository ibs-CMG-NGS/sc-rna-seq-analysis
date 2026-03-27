"""
src/scripts/utils.py
Snakemake 스크립트 공통 유틸리티
"""

import sys
import logging
from pathlib import Path


def get_snakemake_logger(log_path: str) -> logging.Logger:
    """Snakemake log 파일과 stderr 양쪽에 출력하는 logger 반환"""
    logger = logging.getLogger("pipeline")
    logger.setLevel(logging.INFO)

    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")

    # 파일 핸들러
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    # stderr 핸들러
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    return logger


def ensure_output_dirs(*paths: str) -> None:
    """출력 파일 경로의 부모 디렉토리를 생성"""
    for p in paths:
        Path(p).parent.mkdir(parents=True, exist_ok=True)
