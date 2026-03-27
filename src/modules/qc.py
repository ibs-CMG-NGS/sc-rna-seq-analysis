"""
src/modules/qc.py
QC 관련 유틸리티 함수 (notebooks에서 직접 호출용)
"""

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
from typing import Optional


def summarize_qc(adata: ad.AnnData) -> pd.DataFrame:
    """
    QC 지표 요약 통계 반환 (per sample)
    """
    required = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    cols = [c for c in required if c in adata.obs.columns]

    if "sample" in adata.obs.columns:
        return adata.obs.groupby("sample")[cols].describe()
    return adata.obs[cols].describe()


def plot_qc_violin(
    adata: ad.AnnData,
    keys: Optional[list] = None,
    groupby: Optional[str] = "sample",
    save: Optional[str] = None,
) -> None:
    """QC 지표 violin plot"""
    if keys is None:
        keys = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    keys = [k for k in keys if k in adata.obs.columns]

    sc.pl.violin(adata, keys, groupby=groupby, jitter=False, rotation=45, show=False)
    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
    plt.show()


def plot_qc_scatter(
    adata: ad.AnnData,
    save: Optional[str] = None,
) -> None:
    """counts vs genes scatter, mt% scatter"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", ax=axes[0], show=False)
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", ax=axes[1], show=False)
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
    plt.show()
