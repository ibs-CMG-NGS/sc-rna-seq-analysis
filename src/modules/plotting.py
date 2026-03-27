"""
src/modules/plotting.py
재사용 가능한 시각화 함수 (notebooks 공통 사용)
"""

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from typing import Optional, List, Union
import warnings


def umap_panel(
    adata: ad.AnnData,
    color: Union[str, List[str]],
    ncols: int = 4,
    figsize_per_panel: tuple = (4, 3.5),
    save: Optional[str] = None,
    **kwargs,
) -> None:
    """다중 색상 UMAP 패널"""
    if isinstance(color, str):
        color = [color]
    color = [c for c in color if c in adata.obs.columns or c in adata.var_names]

    nrows = int(np.ceil(len(color) / ncols))
    fig, axes = plt.subplots(
        nrows, min(ncols, len(color)),
        figsize=(figsize_per_panel[0] * min(ncols, len(color)), figsize_per_panel[1] * nrows),
    )
    axes = np.array(axes).flatten()

    for i, col in enumerate(color):
        sc.pl.umap(adata, color=col, ax=axes[i], show=False, **kwargs)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
    plt.show()


def marker_dotplot(
    adata: ad.AnnData,
    markers: dict,
    groupby: str = "cell_type",
    save: Optional[str] = None,
    **kwargs,
) -> None:
    """
    마커 유전자 dotplot

    Parameters
    ----------
    markers : dict
        {'cell_type': ['gene1', 'gene2']} 형태
    """
    gene_lists = {}
    for ct, info in markers.items():
        if isinstance(info, dict):
            genes = info.get("markers", [])
        else:
            genes = info
        gene_lists[ct] = [g for g in genes if g in adata.var_names]

    all_genes = [g for genes in gene_lists.values() for g in genes]
    if not all_genes:
        warnings.warn("No marker genes found in adata.var_names")
        return

    sc.pl.dotplot(
        adata,
        var_names=gene_lists,
        groupby=groupby,
        show=False,
        **kwargs,
    )
    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
    plt.show()


def volcano_plot(
    de_df,
    logfc_col: str = "logfoldchanges",
    pval_col: str = "pvals_adj",
    gene_col: str = "names",
    logfc_threshold: float = 1.0,
    pval_threshold: float = 0.05,
    n_label: int = 10,
    save: Optional[str] = None,
) -> None:
    """DEG volcano plot"""
    import pandas as pd
    import numpy as np

    df = de_df.copy()
    df["-log10_pval"] = -np.log10(df[pval_col].clip(1e-300))
    df["significant"] = (abs(df[logfc_col]) >= logfc_threshold) & (df[pval_col] < pval_threshold)
    df["color"] = "grey"
    df.loc[df[logfc_col] >= logfc_threshold, "color"] = "red"
    df.loc[df[logfc_col] <= -logfc_threshold, "color"] = "blue"

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df[logfc_col], df["-log10_pval"], c=df["color"], alpha=0.4, s=8)
    ax.axvline(x=logfc_threshold, color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(x=-logfc_threshold, color="gray", linestyle="--", linewidth=0.8)
    ax.axhline(y=-np.log10(pval_threshold), color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10 Adjusted p-value")

    # 상위 유전자 라벨링
    top = df[df["significant"]].nlargest(n_label, "-log10_pval")
    for _, row in top.iterrows():
        ax.annotate(row[gene_col], (row[logfc_col], row["-log10_pval"]), fontsize=7)

    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
    plt.show()
