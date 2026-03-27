"""
src/modules/annotation.py
Cell type annotation 유틸리티

주요 기능:
  - CellTypist 자동 분류
  - 수동 annotation 적용
  - Marker 기반 점수 계산
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict, List
import warnings


# =============================================================================
# CellTypist 자동 분류
# =============================================================================
def run_celltypist(
    adata: ad.AnnData,
    model: str = "Immune_All_Low.pkl",
    majority_voting: bool = True,
    over_clustering: Optional[str] = "leiden",
) -> ad.AnnData:
    """
    CellTypist으로 자동 cell type 분류

    Parameters
    ----------
    adata : AnnData (log-normalized 상태여야 함)
    model : CellTypist 모델명. 없으면 자동 다운로드.
    majority_voting : True이면 클러스터 내 majority vote로 안정화
    over_clustering : majority_voting에 사용할 obs 칼럼명

    Returns
    -------
    AnnData with 'celltypist_cell_type' and 'celltypist_conf_score' in obs
    """
    try:
        import celltypist
        from celltypist import models
    except ImportError:
        raise ImportError("celltypist not installed. Run: pip install celltypist")

    # 모델 다운로드 (이미 있으면 스킵)
    models.download_models(model=model, force_update=False)

    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=majority_voting,
        over_clustering=over_clustering,
    )
    adata = predictions.to_adata()
    adata.obs["celltypist_cell_type"] = adata.obs["majority_voting"] if majority_voting else adata.obs["predicted_labels"]
    adata.obs["celltypist_conf_score"] = adata.obs["conf_score"]
    return adata


# =============================================================================
# 수동 annotation 적용
# =============================================================================
def apply_manual_annotation(
    adata: ad.AnnData,
    cluster_to_celltype: Dict[str, str],
    cluster_key: str = "leiden",
    output_key: str = "cell_type",
) -> ad.AnnData:
    """
    클러스터 → cell type 수동 매핑 적용

    Parameters
    ----------
    cluster_to_celltype : {'0': 'CD4 T', '1': 'B cell', ...}
    """
    mapping = pd.Series(cluster_to_celltype)
    adata.obs[output_key] = adata.obs[cluster_key].map(mapping)
    n_unmapped = adata.obs[output_key].isna().sum()
    if n_unmapped > 0:
        warnings.warn(f"{n_unmapped} cells have no annotation (NaN). Check cluster_to_celltype mapping.")
    return adata


# =============================================================================
# Marker 기반 점수 계산
# =============================================================================
def score_cell_types(
    adata: ad.AnnData,
    markers: Dict,
    prefix: str = "score_",
) -> ad.AnnData:
    """
    markers.yaml 기반으로 각 cell type 점수 계산 (sc.tl.score_genes)

    Parameters
    ----------
    markers : markers.yaml 로드한 dict
    """
    for category, subtypes in markers.items():
        for subtype, info in subtypes.items():
            if isinstance(info, dict):
                gene_list = info.get("markers", [])
            else:
                gene_list = info

            gene_list = [g for g in gene_list if g in adata.var_names]
            if not gene_list:
                continue

            key = f"{prefix}{subtype}"
            sc.tl.score_genes(adata, gene_list=gene_list, score_name=key)

    return adata
