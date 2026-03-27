# src/rules/hvg.smk
# Highly Variable Gene 선정 및 스케일링 규칙


rule select_hvg:
    """배치를 고려한 HVG 선정 및 스케일링"""
    input:
        h5ad="output/checkpoints/{dataset}/03_normalized.h5ad",
    output:
        h5ad="output/checkpoints/{dataset}/04_hvg.h5ad",
        fig_hvg="output/figures/qc/{dataset}/hvg_plot.png",
    log:
        "output/logs/{dataset}/04_select_hvg.log",
    params:
        flavor=config["hvg"]["flavor"],
        n_top_genes=config["hvg"]["n_top_genes"],
        batch_key=config["hvg"]["batch_key"],
    script:
        "../scripts/select_hvg.py"
