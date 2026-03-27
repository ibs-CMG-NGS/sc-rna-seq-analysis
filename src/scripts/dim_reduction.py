"""
src/scripts/dim_reduction.py
PCA, Harmony batch correction, UMAP 스크립트

Snakemake rules:
  - pca: PCA 계산
  - batch_correct_umap: Harmony → neighbors → UMAP
"""

import sys
sys.path.insert(0, str(__file__).rsplit("/scripts/", 1)[0])

import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from modules.io import load_adata, save_adata
from scripts.utils import get_snakemake_logger, ensure_output_dirs


# =============================================================================
# pca rule
# =============================================================================
def run_pca(snakemake):
    log = get_snakemake_logger(snakemake.log[0])
    ensure_output_dirs(snakemake.output.h5ad, snakemake.output.fig_variance)

    log.info(f"Loading: {snakemake.input.h5ad}")
    adata = load_adata(snakemake.input.h5ad)

    sc.tl.pca(
        adata,
        n_comps=snakemake.params.n_comps,
        random_state=snakemake.params.random_state,
    )
    log.info(f"PCA computed: {snakemake.params.n_comps} components")

    # Elbow plot (explained variance)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
    plt.savefig(snakemake.output.fig_variance, dpi=150, bbox_inches="tight")
    plt.close()

    save_adata(adata, snakemake.output.h5ad)
    log.info(f"Saved: {snakemake.output.h5ad}")


# =============================================================================
# batch_correct_umap rule
# =============================================================================
def run_batch_correct_umap(snakemake):
    log = get_snakemake_logger(snakemake.log[0])
    ensure_output_dirs(snakemake.output.h5ad, snakemake.output.fig_umap_sample)

    log.info(f"Loading: {snakemake.input.h5ad}")
    adata = load_adata(snakemake.input.h5ad)

    harmony_params = snakemake.params.harmony
    neighbors_params = snakemake.params.neighbors
    umap_params = snakemake.params.umap
    batch_key = harmony_params["batch_key"]

    # Harmony batch correction
    if batch_key in adata.obs.columns:
        try:
            import harmonypy as hm
            ho = hm.run_harmony(
                adata.obsm["X_pca"],
                adata.obs,
                batch_key,
                max_iter_harmony=harmony_params.get("max_iter_harmony", 20),
            )
            adata.obsm["X_pca_harmony"] = ho.Z_corr.T
            use_rep = "X_pca_harmony"
            log.info(f"Harmony batch correction applied (batch_key='{batch_key}')")
        except ImportError:
            log.warning("harmonypy not installed. Using raw PCA.")
            use_rep = "X_pca"
    else:
        log.warning(f"batch_key '{batch_key}' not in obs. Skipping Harmony.")
        use_rep = "X_pca"

    # Neighbor graph
    sc.pp.neighbors(
        adata,
        n_neighbors=neighbors_params["n_neighbors"],
        n_pcs=neighbors_params["n_pcs"],
        use_rep=use_rep,
        metric=neighbors_params.get("metric", "euclidean"),
        random_state=umap_params.get("random_state", 42),
    )
    log.info("Neighbor graph computed")

    # UMAP
    sc.tl.umap(
        adata,
        min_dist=umap_params["min_dist"],
        spread=umap_params["spread"],
        random_state=umap_params.get("random_state", 42),
    )
    log.info("UMAP computed")

    # UMAP colored by sample
    color_key = batch_key if batch_key in adata.obs.columns else None
    if color_key:
        sc.pl.umap(adata, color=color_key, show=False)
        plt.savefig(snakemake.output.fig_umap_sample, dpi=150, bbox_inches="tight")
        plt.close()
    else:
        # 빈 figure 저장
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No sample metadata", ha="center", va="center")
        plt.savefig(snakemake.output.fig_umap_sample, dpi=150)
        plt.close()

    save_adata(adata, snakemake.output.h5ad)
    log.info(f"Saved: {snakemake.output.h5ad}")


# =============================================================================
# Snakemake 진입점
# =============================================================================
if "snakemake" in dir():
    rule_name = snakemake.rule  # type: ignore[name-defined]
    if rule_name == "pca":
        run_pca(snakemake)  # type: ignore[name-defined]
    elif rule_name == "batch_correct_umap":
        run_batch_correct_umap(snakemake)  # type: ignore[name-defined]
    else:
        raise ValueError(f"Unknown rule: {rule_name}")
