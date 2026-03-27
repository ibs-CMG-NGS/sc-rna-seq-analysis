# src/rules/normalize.smk
# 정규화 및 log1p 변환 규칙


rule normalize:
    """총합 정규화 + log1p 변환. 원본 count는 adata.raw에 보존."""
    input:
        h5ad="output/checkpoints/{dataset}/02_filtered.h5ad",
    output:
        h5ad="output/checkpoints/{dataset}/03_normalized.h5ad",
    log:
        "output/logs/{dataset}/03_normalize.log",
    params:
        target_sum=config["normalization"]["target_sum"],
    script:
        "../scripts/normalize.py"
