"""
src/scripts/cluster.py
Leiden 클러스터링 스크립트

Snakemake rule: cluster
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
    ensure_output_dirs(  # type: ignore[name-defined]
        snakemake.output.h5ad,
        snakemake.output.fig_umap_clusters,
        snakemake.output.fig_umap_qc,
    )

    log.info(f"Loading: {snakemake.input.h5ad}")  # type: ignore[name-defined]
    adata = load_adata(snakemake.input.h5ad)  # type: ignore[name-defined]

    p = snakemake.params  # type: ignore[name-defined]

    # 다중 resolution으로 클러스터링
    for res in p.resolutions:
        key = f"leiden_{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key, random_state=p.random_state)
        n_clusters = adata.obs[key].nunique()
        log.info(f"  resolution={res}: {n_clusters} clusters")

    # 기본 resolution을 'leiden'으로 설정
    default_key = f"leiden_{p.default_resolution}"
    adata.obs["leiden"] = adata.obs[default_key].copy()
    log.info(f"Default clustering: 'leiden' = leiden_{p.default_resolution}")

    # Dendrogram
    sc.tl.dendrogram(adata, groupby="leiden")

    # UMAP colored by clusters
    sc.pl.umap(adata, color="leiden", legend_loc="on data", show=False)
    plt.savefig(snakemake.output.fig_umap_clusters, dpi=150, bbox_inches="tight")
    plt.close()

    # UMAP colored by QC metrics
    qc_colors = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    qc_colors = [c for c in qc_colors if c in adata.obs.columns]
    if qc_colors:
        sc.pl.umap(adata, color=qc_colors, ncols=len(qc_colors), show=False)
        plt.savefig(snakemake.output.fig_umap_qc, dpi=150, bbox_inches="tight")
        plt.close()
    else:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No QC metrics found", ha="center", va="center")
        plt.savefig(snakemake.output.fig_umap_qc, dpi=150)
        plt.close()

    save_adata(adata, snakemake.output.h5ad)  # type: ignore[name-defined]
    log.info(f"Saved: {snakemake.output.h5ad}")  # type: ignore[name-defined]
