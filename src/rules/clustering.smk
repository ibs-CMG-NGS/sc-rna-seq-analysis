# src/rules/clustering.smk
# Leiden 클러스터링 규칙


rule cluster:
    """다중 resolution으로 Leiden 클러스터링 수행"""
    input:
        h5ad="output/checkpoints/{dataset}/06_umap.h5ad",
    output:
        h5ad="output/checkpoints/{dataset}/07_clustered.h5ad",
        fig_umap_clusters="output/figures/umap/{dataset}/umap_clusters.png",
        fig_umap_qc="output/figures/umap/{dataset}/umap_qc_metrics.png",
    log:
        "output/logs/{dataset}/07_cluster.log",
    params:
        algorithm=config["clustering"]["algorithm"],
        resolutions=config["clustering"]["resolutions"],
        default_resolution=config["clustering"]["default_resolution"],
        random_state=config["clustering"]["random_state"],
    script:
        "../scripts/cluster.py"
