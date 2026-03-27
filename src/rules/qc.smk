# src/rules/qc.smk
# QC 지표 계산 및 세포 필터링 규칙


rule compute_qc:
    """원본 h5ad에 QC 지표(mt%, rb%, n_genes, n_counts)를 계산하고 저장"""
    input:
        h5ad=lambda wc: config["datasets"][wc.dataset]["path"],
    output:
        h5ad="output/checkpoints/{dataset}/01_qc.h5ad",
        fig_violin="output/figures/qc/{dataset}/violin_raw.png",
        fig_scatter="output/figures/qc/{dataset}/scatter_qc.png",
    log:
        "output/logs/{dataset}/01_compute_qc.log",
    params:
        mt_prefix=config["genome"]["mt_prefix"],
        rb_prefixes=config["genome"]["rb_prefixes"],
        fmt=lambda wc: config["datasets"][wc.dataset].get("format", "auto"),
    script:
        "../scripts/qc_filter.py"


rule filter_cells:
    """QC 임계값으로 세포 필터링 및 Scrublet doublet 제거"""
    input:
        h5ad="output/checkpoints/{dataset}/01_qc.h5ad",
    output:
        h5ad="output/checkpoints/{dataset}/02_filtered.h5ad",
        fig_violin="output/figures/qc/{dataset}/violin_filtered.png",
        stats="output/tables/cell_counts/{dataset}_filter_stats.tsv",
    log:
        "output/logs/{dataset}/02_filter_cells.log",
    params:
        min_genes=config["qc"]["min_genes"],
        max_genes=config["qc"]["max_genes"],
        min_counts=config["qc"]["min_counts"],
        max_pct_mt=config["qc"]["max_pct_mt"],
        scrublet=config["qc"]["scrublet"],
    script:
        "../scripts/qc_filter.py"
