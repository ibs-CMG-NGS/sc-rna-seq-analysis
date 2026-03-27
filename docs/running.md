# Pipeline Running Guide

## 기본 실행

```bash
conda activate sc-rna-seq-analysis

# 프로젝트 루트에서 실행
cd /home/ygkim/sc-rna-seq-analysis

# 실행 계획 확인 (실제 실행 없음)
snakemake --snakefile src/Snakefile --dryrun

# 실행 (4 cores)
snakemake --snakefile src/Snakefile --cores 4
```

---

## 자주 쓰는 Snakemake 옵션

| 옵션 | 설명 |
|------|------|
| `--dryrun` | 실행하지 않고 어떤 규칙이 실행될지 확인 |
| `--cores 4` | 병렬 실행 core 수 |
| `--forcerun rule_name` | 특정 규칙 강제 재실행 |
| `--until rule_name` | 특정 규칙까지만 실행 |
| `--rerun-incomplete` | 불완전하게 종료된 작업 재실행 |
| `--keep-going` | 일부 규칙 실패해도 다른 규칙 계속 실행 |
| `--report report.html` | 실행 결과 HTML 보고서 생성 |

---

## 특정 데이터셋만 실행

```bash
# human_genes_only만 전체 파이프라인 실행
snakemake --snakefile src/Snakefile --cores 4 \
  output/checkpoints/human_genes_only/07_clustered.h5ad

# EGFP 데이터셋 QC까지만
snakemake --snakefile src/Snakefile --cores 4 \
  output/checkpoints/egfp_genes_only/02_filtered.h5ad
```

---

## 파라미터 변경 후 재실행

`configs/pipeline.yaml`에서 파라미터 변경 시 영향받는 단계만 자동 재실행됩니다.

```bash
# 예: min_genes를 200 → 300으로 변경 후
# filter_cells 이후 모든 단계가 자동 재실행됨
snakemake --snakefile src/Snakefile --cores 4
```

---

## 특정 단계 강제 재실행

```bash
# normalize 단계부터 다시 실행
snakemake --snakefile src/Snakefile --cores 4 \
  --forcerun normalize
```

---

## Jupyter Notebook 실행

Snakemake 파이프라인 완료 후 notebook 순서대로 실행:

```bash
conda activate sc-rna-seq-analysis
cd /home/ygkim/sc-rna-seq-analysis

# JupyterLab 실행
jupyter lab --notebook-dir=src/notebooks

# 또는 VSCode에서 notebook 파일 직접 열기
```

Notebook 실행 순서:
1. `01_cluster_review.ipynb`
2. `02_cell_type_annotation.ipynb`
3. `03_deg_analysis.ipynb`
4. `04_trajectory_paga.ipynb`
5. `05_figures.ipynb`

---

## 로그 확인

```bash
# 특정 단계 로그
cat output/logs/human_genes_only/03_normalize.log

# 모든 로그 실시간 확인
tail -f output/logs/human_genes_only/*.log
```

---

## 새 데이터셋 추가

1. `configs/samples.yaml`에 데이터셋 추가 (`run_in_pipeline: true`)
2. `snakemake --snakefile src/Snakefile --dryrun` 으로 확인
3. `snakemake --snakefile src/Snakefile --cores 4` 실행
