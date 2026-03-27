"""
src/scripts/qc_filter.py
QC 지표 계산 및 세포 필터링 스크립트

Snakemake rules:
  - compute_qc: QC 지표 추가 (rule name: compute_qc)
  - filter_cells: 임계값 기반 필터링 + Scrublet (rule name: filter_cells)
"""

import sys
sys.path.insert(0, str(__file__).rsplit("/scripts/", 1)[0])  # src/ 경로 추가

import scanpy as sc
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from modules.io import load_adata, save_adata
from scripts.utils import get_snakemake_logger, ensure_output_dirs


# =============================================================================
# compute_qc rule
# =============================================================================
def run_compute_qc(snakemake):
    log = get_snakemake_logger(snakemake.log[0])
    ensure_output_dirs(snakemake.output.h5ad, snakemake.output.fig_violin, snakemake.output.fig_scatter)

    fmt = snakemake.params.get("fmt", "auto")
    log.info(f"Loading: {snakemake.input.h5ad} (format={fmt})")
    adata = load_adata(snakemake.input.h5ad, fmt=fmt)
    log.info(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    mt_prefix = snakemake.params.mt_prefix
    rb_prefixes = snakemake.params.rb_prefixes

    # 미토콘드리아 / 리보솜 유전자 플래그
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    adata.var["rb"] = adata.var_names.str.startswith(tuple(rb_prefixes))

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "rb"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    log.info("QC metrics calculated")

    # 기존 split-pipe 클러스터 보존
    if "leiden" in adata.obs.columns:
        adata.obs["leiden_splitpipe"] = adata.obs["leiden"].copy()
        log.info("Preserved original split-pipe leiden clusters as 'leiden_splitpipe'")

    # Violin plot
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    plt.savefig(snakemake.output.fig_violin, dpi=150, bbox_inches="tight")
    plt.close()

    # Scatter plot
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", ax=axes[0], show=False)
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", ax=axes[1], show=False)
    plt.tight_layout()
    plt.savefig(snakemake.output.fig_scatter, dpi=150, bbox_inches="tight")
    plt.close()

    save_adata(adata, snakemake.output.h5ad)
    log.info(f"Saved: {snakemake.output.h5ad}")


# =============================================================================
# filter_cells rule
# =============================================================================
def run_filter_cells(snakemake):
    log = get_snakemake_logger(snakemake.log[0])
    ensure_output_dirs(snakemake.output.h5ad, snakemake.output.fig_violin, snakemake.output.stats)

    log.info(f"Loading: {snakemake.input.h5ad}")
    adata = load_adata(snakemake.input.h5ad)
    n_before = adata.n_obs
    log.info(f"Before filtering: {n_before} cells")

    p = snakemake.params

    # 기본 임계값 필터
    sc.pp.filter_cells(adata, min_genes=p.min_genes)
    sc.pp.filter_cells(adata, max_genes=p.max_genes)
    sc.pp.filter_cells(adata, min_counts=p.min_counts)
    adata = adata[adata.obs.pct_counts_mt < p.max_pct_mt].copy()
    log.info(f"After threshold filtering: {adata.n_obs} cells")

    # Scrublet doublet detection
    try:
        import scrublet as scr
        scrub = scr.Scrublet(
            adata.X,
            expected_doublet_rate=p.scrublet["expected_doublet_rate"],
        )
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
        adata.obs["doublet_score"] = doublet_scores
        adata.obs["predicted_doublet"] = predicted_doublets

        # 임계값 오버라이드 (설정값 사용)
        threshold = p.scrublet["doublet_score_threshold"]
        adata.obs["predicted_doublet"] = adata.obs["doublet_score"] > threshold

        n_doublets = adata.obs["predicted_doublet"].sum()
        log.info(f"Scrublet detected {n_doublets} doublets (threshold={threshold})")
        adata = adata[~adata.obs["predicted_doublet"]].copy()
        log.info(f"After doublet removal: {adata.n_obs} cells")
    except Exception as e:
        log.warning(f"Scrublet failed: {e}. Skipping doublet removal.")

    # 유전자 필터 (최소 3개 세포에서 발현)
    sc.pp.filter_genes(adata, min_cells=3)
    log.info(f"After gene filtering: {adata.n_vars} genes")

    # 필터링 통계 저장
    stats = pd.DataFrame({
        "metric": ["cells_before", "cells_after", "genes_after"],
        "value": [n_before, adata.n_obs, adata.n_vars],
    })
    stats.to_csv(snakemake.output.stats, sep="\t", index=False)

    # Violin plot (필터링 후)
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    plt.savefig(snakemake.output.fig_violin, dpi=150, bbox_inches="tight")
    plt.close()

    save_adata(adata, snakemake.output.h5ad)
    log.info(f"Saved: {snakemake.output.h5ad}")


# =============================================================================
# Snakemake 진입점 — rule 이름으로 분기
# =============================================================================
if "snakemake" in dir():
    rule_name = snakemake.rule  # type: ignore[name-defined]
    if rule_name == "compute_qc":
        run_compute_qc(snakemake)  # type: ignore[name-defined]
    elif rule_name == "filter_cells":
        run_filter_cells(snakemake)  # type: ignore[name-defined]
    else:
        raise ValueError(f"Unknown rule: {rule_name}")
