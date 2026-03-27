# scRNA-seq Tertiary Analysis Pipeline

Single-cell RNA-seq downstream analysis pipeline using **Scanpy + Snakemake**.

## Overview

| Item | Detail |
|------|--------|
| Input | AnnData h5ad (split-pipe), 10x MTX, 10x H5, loom |
| Main library | [Scanpy](https://scanpy.readthedocs.io) |
| Workflow manager | Snakemake |
| Batch correction | Harmony |
| Cell annotation | CellTypist |
| Trajectory | PAGA + DPT + CytoTRACE / scVelo (auto-detect) |
| Environment | conda `sc-rna-seq-analysis` (Python 3.10) |

---

## Directory Structure

```
sc-rna-seq-analysis/
├── environment.yml          # Conda 환경 정의
├── configs/
│   ├── pipeline.yaml        # QC, 정규화, 클러스터링 파라미터 (template)
│   ├── samples.yaml         # 입력 데이터셋 및 샘플 메타데이터 (template)
│   └── markers.yaml         # Cell type 마커 유전자 목록
├── docs/
│   ├── setup.md             # 설치 가이드
│   ├── configuration.md     # 파라미터 레퍼런스
│   ├── running.md           # 실행 방법
│   ├── notebooks_guide.md   # Notebook 사용 가이드
│   └── roadmap.md           # 향후 기능 계획
├── src/
│   ├── Snakefile            # 메인 워크플로우
│   ├── rules/               # Snakemake rule 파일 (단계별)
│   ├── scripts/             # Snakemake가 호출하는 Python 스크립트
│   ├── modules/             # Notebooks/scripts 공통 Python 라이브러리
│   └── notebooks/           # Jupyter notebooks (탐색적 분석)
├── data/                    # 입력 데이터 (gitignore, 구조만 유지)
└── output/                  # 분석 결과 (gitignore, 구조만 유지)
```

---

## Quick Start

### 1. 환경 설치

```bash
conda env create -f environment.yml
conda activate sc-rna-seq-analysis

# Jupyter kernel 등록
python -m ipykernel install --user --name sc-rna-seq-analysis --display-name "sc-rna-seq-analysis"
```

### 2. 데이터 및 설정 준비

```bash
# 데이터 배치
cp /path/to/your/data.h5ad data/

# config template 복사 후 수정 (gitignore 대상)
cp configs/samples.yaml configs/_samples.yaml
# configs/_samples.yaml 에서 path, condition 등 수정
```

> `configs/_*.yaml` 파일은 gitignore 대상입니다. 민감한 경로/조건 정보를 안전하게 관리할 수 있습니다.

### 3. 파이프라인 실행

```bash
# 실행 계획 확인 (dry-run)
snakemake --snakefile src/Snakefile --dryrun

# 실행
snakemake --snakefile src/Snakefile --cores 4
```

### 4. Jupyter Notebooks 실행

Snakemake 완료 후 순서대로 실행:

```bash
jupyter lab --notebook-dir=src/notebooks
```

| Notebook | 목적 |
|----------|------|
| `01_cluster_review.ipynb` | QC 검토, 저품질 클러스터 제거 |
| `02_cell_type_annotation.ipynb` | CellTypist 자동 분류 + 수동 보정 |
| `03_deg_analysis.ipynb` | DEG 분석 (Wilcoxon), volcano plot |
| `04_trajectory_paga.ipynb` | Trajectory (PAGA + DPT / scVelo) |
| `05_figures.ipynb` | 출판용 figure 생성 |

---

## Pipeline Steps

```
Raw h5ad
   │
   ├─ [1] compute_qc       mt%, rb%, n_genes, n_counts 계산
   ├─ [2] filter_cells      임계값 필터링 + Scrublet doublet 제거
   ├─ [3] normalize         normalize_total + log1p (raw 보존)
   ├─ [4] select_hvg        HVG 선정 (seurat_v3) + 스케일링
   ├─ [5] pca               PCA (50 components)
   ├─ [6] batch_correct_umap  Harmony batch correction + UMAP
   └─ [7] cluster           Leiden 클러스터링 (다중 resolution)
```

---

## Trajectory Analysis

BAM/loom 파일 유무에 따라 자동 분기:

```
loom 또는 spliced layer 존재  →  scVelo (RNA velocity)
BAM 파일 존재                 →  velocyto 실행 안내 → scVelo
없음                          →  PAGA + DPT + CytoTRACE
```

scVelo 사용 시 추가 설치:
```bash
pip install scvelo velocyto
```

---

## Supported Input Formats

| Format | Source | Path example |
|--------|--------|-------------|
| `h5ad` | split-pipe, alevin | `data/sample.h5ad` |
| `10x_mtx` | CellRanger, STARsolo | `data/sample/` (directory) |
| `10x_h5` | CellRanger | `data/filtered_feature_bc_matrix.h5` |
| `loom` | velocyto | `data/sample.loom` |

`configs/samples.yaml`에서 `format: auto`로 설정하면 자동 감지합니다.

---

## Adding New Analysis Modules

```bash
# 새 Snakemake 단계 추가
touch src/rules/new_step.smk
# src/Snakefile에 include: "rules/new_step.smk" 추가

# 새 탐색 분석 추가
touch src/notebooks/06_new_analysis.ipynb

# 새 데이터셋 추가
# configs/samples.yaml 에 항목 추가 후 snakemake 실행
```

---

## Documentation

- [Setup Guide](docs/setup.md)
- [Configuration Reference](docs/configuration.md)
- [Running the Pipeline](docs/running.md)
- [Notebooks Guide](docs/notebooks_guide.md)
- [Roadmap](docs/roadmap.md)
