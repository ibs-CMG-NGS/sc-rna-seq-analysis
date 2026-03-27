# src/rules/dim_reduction.smk
# PCA → Harmony batch correction → UMAP 규칙


rule pca:
    """PCA 계산 (HVG subset 사용)"""
    input:
        h5ad="output/checkpoints/{dataset}/04_hvg.h5ad",
    output:
        h5ad="output/checkpoints/{dataset}/05_pca.h5ad",
        fig_variance="output/figures/qc/{dataset}/pca_variance.png",
    log:
        "output/logs/{dataset}/05_pca.log",
    params:
        n_comps=config["pca"]["n_comps"],
        random_state=config["umap"]["random_state"],
    script:
        "../scripts/dim_reduction.py"


rule batch_correct_umap:
    """Harmony batch correction 후 neighbor graph 및 UMAP 계산"""
    input:
        h5ad="output/checkpoints/{dataset}/05_pca.h5ad",
    output:
        h5ad="output/checkpoints/{dataset}/06_umap.h5ad",
        fig_umap_sample="output/figures/umap/{dataset}/umap_by_sample.png",
    log:
        "output/logs/{dataset}/06_umap.log",
    params:
        harmony=config["harmony"],
        neighbors=config["neighbors"],
        umap=config["umap"],
    script:
        "../scripts/dim_reduction.py"
