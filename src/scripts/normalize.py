"""
src/scripts/normalize.py
정규화 및 log1p 변환 스크립트

Snakemake rule: normalize
"""

import sys
sys.path.insert(0, str(__file__).rsplit("/scripts/", 1)[0])

import scanpy as sc

from modules.io import load_adata, save_adata
from scripts.utils import get_snakemake_logger, ensure_output_dirs


if "snakemake" in dir():
    log = get_snakemake_logger(snakemake.log[0])  # type: ignore[name-defined]
    ensure_output_dirs(snakemake.output.h5ad)  # type: ignore[name-defined]

    log.info(f"Loading: {snakemake.input.h5ad}")  # type: ignore[name-defined]
    adata = load_adata(snakemake.input.h5ad)  # type: ignore[name-defined]

    # 원본 count 보존 (downstream DEG, velocity 등에 필요)
    adata.raw = adata

    target_sum = snakemake.params.target_sum  # type: ignore[name-defined]
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    log.info(f"Normalized to {target_sum} counts/cell and log1p transformed")

    save_adata(adata, snakemake.output.h5ad)  # type: ignore[name-defined]
    log.info(f"Saved: {snakemake.output.h5ad}")  # type: ignore[name-defined]
