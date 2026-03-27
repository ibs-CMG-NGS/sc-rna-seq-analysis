"""
src/scripts/select_hvg.py
Highly Variable Gene 선정 및 스케일링 스크립트

Snakemake rule: select_hvg
"""

import sys
sys.path.insert(0, str(__file__).rsplit("/scripts/", 1)[0])

import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from modules.io import load_adata, save_adata
from scripts.utils import get_snakemake_logger, ensure_output_dirs


if "snakemake" in dir():
    log = get_snakemake_logger(snakemake.log[0])  # type: ignore[name-defined]
    ensure_output_dirs(snakemake.output.h5ad, snakemake.output.fig_hvg)  # type: ignore[name-defined]

    log.info(f"Loading: {snakemake.input.h5ad}")  # type: ignore[name-defined]
    adata = load_adata(snakemake.input.h5ad)  # type: ignore[name-defined]

    p = snakemake.params  # type: ignore[name-defined]

    # batch_key가 obs에 있는지 확인
    batch_key = p.batch_key if p.batch_key in adata.obs.columns else None
    if batch_key is None:
        log.warning(f"batch_key '{p.batch_key}' not found in obs. Running HVG without batch correction.")

    sc.pp.highly_variable_genes(
        adata,
        flavor=p.flavor,
        n_top_genes=p.n_top_genes,
        batch_key=batch_key,
    )
    n_hvg = adata.var["highly_variable"].sum()
    log.info(f"Selected {n_hvg} highly variable genes")

    # HVG dispersion plot
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(snakemake.output.fig_hvg, dpi=150, bbox_inches="tight")
    plt.close()

    # HVG subset으로 제한 후 스케일링
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    log.info(f"Scaled data (max_value=10). Shape: {adata.shape}")

    save_adata(adata, snakemake.output.h5ad)  # type: ignore[name-defined]
    log.info(f"Saved: {snakemake.output.h5ad}")  # type: ignore[name-defined]
