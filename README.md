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
│   ├── pipeline.yaml        # QC, 정규화, 클러스터링, velocity 파라미터
│   ├── samples.yaml         # 입력 데이터셋, 샘플 메타데이터, BAM/loom 경로
│   └── markers.yaml         # Cell type 마커 유전자 목록
├── docs/
│   ├── setup.md             # 설치 가이드
│   ├── configuration.md     # 파라미터 레퍼런스
│   ├── running.md           # 실행 방법
│   ├── notebooks_guide.md   # Notebook 사용 가이드
│   └── roadmap.md           # 향후 기능 계획
├── src/
│   ├── Snakefile            # 메인 워크플로우
│   ├── rules/
│   │   ├── qc.smk           # QC 및 필터링
│   │   ├── normalize.smk    # 정규화
│   │   ├── hvg.smk          # HVG 선정
│   │   ├── dim_reduction.smk # PCA, UMAP, Harmony
│   │   ├── clustering.smk   # Leiden 클러스터링
│   │   └── velocity.smk     # velocyto (BAM → loom) — 서버 실행
│   ├── scripts/             # Snakemake가 호출하는 Python 스크립트
│   ├── modules/             # Notebooks/scripts 공통 Python 라이브러리
│   └── notebooks/           # Jupyter notebooks (탐색적 분석)
├── data/                    # 입력 데이터 (gitignore, 구조만 유지)
└── output/
    ├── checkpoints/         # 단계별 h5ad 체크포인트
    ├── figures/             # QC, UMAP, trajectory 시각화
    └── velocity/
        ├── bam/             # 정렬된 BAM (samtools 출력)
        ├── barcodes/        # 샘플별 바코드 TSV (velocyto 입력)
        └── loom/            # velocyto 출력 loom 파일
```

---

## 실행 환경 분리 (로컬 vs 서버)

> 파이프라인 전체를 서버에서 돌릴 수도 있지만, 현재 데이터 배치 기준으로 아래와 같이 분리됩니다.

### 서버에서만 가능한 작업

| 단계 | 이유 |
|------|------|
| `velocyto run` (BAM → loom) | GTF 파일이 서버에 있음 + BAM 파일 대용량 처리 |
| `samtools sort/index` | BAM 전처리 (velocyto 전 필수) |

```bash
# 서버에서 실행
samtools sort -@ 4 -o barcode_headAligned_anno.sorted.bam \
    data/all-sample/barcode_headAligned_anno.bam
samtools index barcode_headAligned_anno.sorted.bam

velocyto run \
    -b all_barcodes.tsv \
    -m /path/to/hg38_rmsk.gtf \        # 선택사항 (권장)
    -o output/velocity/loom/ \
    barcode_headAligned_anno.sorted.bam \
    /path/to/genes.gtf

# 완료 후 loom 파일만 로컬로 복사
scp server:/path/to/output/velocity/loom/*.loom data/loom/
```

### 로컬에서 가능한 작업

| 단계 | 파일 | 비고 |
|------|------|------|
| Snakemake 전처리 (QC → 클러스터링) | h5ad (이미 로컬) | 데이터 크기에 따라 서버 권장 |
| Notebook 01–03 | h5ad | QC 검토, 어노테이션, DEG |
| Notebook 04 (scVelo) | loom + h5ad | loom을 서버에서 복사해온 뒤 실행 |
| Notebook 05 (figures) | h5ad | 출판용 figure 생성 |

> `scVelo dynamical` 모드(`recover_dynamics`)는 수만 세포 기준 수 시간 소요.
> 서버에서 실행 후 결과 h5ad만 로컬로 가져오는 방식도 가능합니다.

### 전체 파이프라인을 서버에서 돌리는 경우

```bash
# 서버에서 전체 실행
snakemake --snakefile src/Snakefile --cores 8

# velocyto까지 포함 (configs/pipeline.yaml에 gtf 경로 설정 후)
snakemake --snakefile src/Snakefile --cores 8 \
    output/velocity/loom/human_genes_only/merged.loom

# Jupyter는 서버에서 실행 후 포트포워딩으로 접속
jupyter lab --no-browser --port=8888
# 로컬에서: ssh -L 8888:localhost:8888 user@server
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
| `04_trajectory_paga.ipynb` | PAGA + DPT + CytoTRACE / **scVelo 완전 분석** (아래 참조) |
| `05_figures.ipynb` | 출판용 figure 생성 |

---

## Pipeline Steps

```
Raw h5ad
   │
   ├─ [1] compute_qc          mt%, rb%, n_genes, n_counts 계산
   ├─ [2] filter_cells         임계값 필터링 + Scrublet doublet 제거
   ├─ [3] normalize            normalize_total + log1p (raw 보존)
   ├─ [4] select_hvg           HVG 선정 (seurat_v3) + 스케일링
   ├─ [5] pca                  PCA (50 components)
   ├─ [6] batch_correct_umap   Harmony batch correction + UMAP
   └─ [7] cluster              Leiden 클러스터링 (다중 resolution)

BAM 파일 (서버) ── velocity.smk ──────────────────────────────────
   │
   ├─ [V1] bam_sort_index      samtools sort + index
   ├─ [V2] extract_barcodes    h5ad → 샘플별 바코드 TSV 추출
   ├─ [V3] velocyto_run        BAM + barcodes + GTF → per-sample loom
   └─ [V4] merge_loom          샘플별 loom → merged.loom (scVelo 입력)
```

---

## Trajectory Analysis (`04_trajectory_paga.ipynb`)

BAM/loom 파일 유무에 따라 자동 분기:

```
loom 또는 spliced layer 존재  →  scVelo 완전 분석 (아래)
BAM 파일 존재                 →  velocyto Snakemake 룰 실행 → scVelo
없음                          →  PAGA + DPT + CytoTRACE
```

### scVelo 분석 내용 (loom 파일 준비 후)

| 단계 | 출력 |
|------|------|
| loom 로드 + adata 병합 | spliced/unspliced layer 확인 |
| `filter_and_normalize` + `moments` | scVelo 전처리 |
| `recover_dynamics` + `velocity` | dynamical 모델 (stochastic 선택 가능) |
| **Stream plot** | 흐름 방향 화살표 UMAP (주 시각화) |
| **Arrow plot** | 개별 세포 화살표 UMAP |
| **Grid plot** | 격자 화살표 UMAP |
| **Latent time** | RNA velocity 기반 pseudotime |
| **Velocity confidence/speed** | velocity 신뢰도 확인 |
| **Phase portrait** | driver gene spliced vs unspliced |

> velocity 파라미터(`scvelo_mode`, `n_top_genes`, `n_jobs` 등)는 `configs/pipeline.yaml`의 `velocity` 섹션에서 설정합니다.

scVelo 사용 시 추가 설치:
```bash
pip install scvelo velocyto loompy
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

- [scRNA-seq 분석 생태계 정리](docs/scrna_seq_ecosystem.md) — 전체 생태계 개념 정리 (플랫폼, 도구, 포맷, 변환)
- [Setup Guide](docs/setup.md)
- [Configuration Reference](docs/configuration.md)
- [Running the Pipeline](docs/running.md)
- [Notebooks Guide](docs/notebooks_guide.md)
- [Roadmap](docs/roadmap.md)
