# scRNA-seq 분석 생태계 정리

## 1. 전체 파이프라인 구조

```
┌─────────────────────────────────────────────────────────────────┐
│                     Layer 1: 실험 (Wet Lab)                      │
│                                                                 │
│   세포 분리 → cDNA 라이브러리 제작 → 시퀀싱 (NGS)                  │
│                                                                 │
│   플랫폼:  10x Genomics    Parse Biosciences    Drop-seq 등      │
└──────────────────────────┬──────────────────────────────────────┘
                           │ FASTQ 파일
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                  Layer 2: 1차 분석 (Primary Analysis)             │
│          FASTQ → 정렬 → 바코드 디멀티플렉싱 → Count Matrix          │
│                                                                 │
│  10x 플랫폼   →  CellRanger   →  MTX 디렉토리 또는 .h5            │
│  Parse 플랫폼 →  split-pipe   →  MTX 또는 .h5ad                  │
│  플랫폼 무관  →  STARsolo     →  MTX 디렉토리                     │
│  플랫폼 무관  →  alevin       →  .h5ad                           │
│  플랫폼 무관  →  kallisto/kb  →  .h5ad 또는 loom                 │
└──────────────────────────┬──────────────────────────────────────┘
                           │ Count Matrix (세포 × 유전자)
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│              Layer 3: 데이터 컨테이너 (Data Container)             │
│                  "분석 결과를 담는 그릇"                            │
│                                                                 │
│     Python 생태계                   R 생태계                      │
│   ┌─────────────────┐             ┌─────────────────┐           │
│   │    AnnData      │  ←── 변환 ──▶│  Seurat Object  │           │
│   │  (저장: .h5ad)  │  SeuratDisk  │  (저장: .rds)   │           │
│   └─────────────────┘ zellkonverter└─────────────────┘           │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│           Layer 4: 2~3차 분석 프레임워크 (Analysis Tools)           │
│                                                                 │
│     Python 생태계                   R 생태계                      │
│   ┌─────────────────┐             ┌─────────────────┐           │
│   │     Scanpy      │             │     Seurat      │           │
│   │  (AnnData 사용) │             │ (Seurat Obj 사용)│           │
│   └─────────────────┘             └─────────────────┘           │
│                                                                 │
│   QC → 정규화 → HVG → PCA → UMAP → 클러스터링 → 어노테이션 → DEG  │
└─────────────────────────────────────────────────────────────────┘
```

---

## 2. 실험 플랫폼

### 주요 플랫폼 비교

| 플랫폼 | 회사 | 방식 | 특징 |
|--------|------|------|------|
| **10x Genomics** | 10x Genomics | Droplet (GEM) | 가장 보편적, 세포당 비용 낮음, 대규모 적용 용이 |
| **Parse Biosciences** | Parse Biosciences | Combinatorial indexing (SPLiT-seq 기반) | 고통량, 고정 세포 사용 가능, 플레이트 기반 |
| **Drop-seq** | 학계 (Harvard) | Droplet | 오픈소스, 자체 구축 가능 |
| **inDrops** | 학계 | Droplet | 오픈소스 |
| **Smart-seq2/3** | 학계 | Full-length | 전장 전사체 커버, 낮은 통량 |

### 플랫폼 선택 기준

```
세포 수가 많고 (>10,000) 비용 절감 필요  →  10x Genomics
고정 조직 또는 핵 분리 샘플              →  Parse Biosciences
전장 전사체 정보 필요 (splicing 분석 등) →  Smart-seq2/3
자체 구축 및 커스터마이징 필요            →  Drop-seq
```

---

## 3. 1차 분석 도구

### 도구별 특징

| 도구 | 연결 플랫폼 | 주요 출력 포맷 | 비고 |
|------|-----------|--------------|------|
| **CellRanger** | 10x 전용 | MTX 디렉토리, `.h5` | 10x 공식 도구, 가장 안정적 |
| **split-pipe** | Parse 전용 | MTX, `.h5ad` | Parse 공식 도구 |
| **STARsolo** | 플랫폼 무관 | MTX 디렉토리 | STAR 기반, 빠르고 유연 |
| **alevin / alevin-fry** | 플랫폼 무관 | `.h5ad` | salmon 기반, 빠른 처리 |
| **kallisto / bustools (kb)** | 플랫폼 무관 | `.h5ad`, loom | 경량, 빠른 처리 |

### 1차 분석 출력물 구조

**MTX 디렉토리 (10x 표준 포맷)**
```
filtered_feature_bc_matrix/
├── barcodes.tsv.gz    ← 세포 바코드 목록
├── features.tsv.gz    ← 유전자 목록 (gene_id, gene_name)
└── matrix.mtx.gz      ← 희소 행렬 (Sparse Matrix)
```

**H5 파일 (CellRanger HDF5)**
```
filtered_feature_bc_matrix.h5   ← 단일 파일, 위 내용 압축
```

**H5AD 파일 (AnnData HDF5)**
```
sample.h5ad   ← AnnData 객체 저장 형식 (Python 생태계 표준)
```

---

## 4. 데이터 컨테이너

### AnnData (Python 생태계)

scRNA-seq 데이터를 구조화된 방식으로 담는 Python 클래스.
`.h5ad` 형식으로 디스크에 저장.

```
AnnData (n_obs × n_vars)
│
├── .X           [cells × genes] count matrix (희소 행렬)
├── .obs         [cells × metadata] 세포 메타데이터
│                 예: sample, cluster, cell_type, n_genes, pct_mt
├── .var         [genes × metadata] 유전자 메타데이터
│                 예: highly_variable, mt, mean_counts
├── .obsm        세포 임베딩 (딕셔너리)
│                 예: X_pca, X_umap, X_harmony
├── .obsp        세포 간 그래프 (희소 행렬)
│                 예: connectivities, distances
├── .uns         비정형 정보 (딕셔너리)
│                 예: leiden_colors, rank_genes_groups, paga
├── .layers      추가 행렬 레이어 (딕셔너리)
│                 예: counts (raw), spliced, unspliced, normalized
└── .raw         정규화 전 원본 count 보존 (adata.raw)
```

### Seurat Object (R 생태계)

R에서 사용하는 동일 목적의 데이터 컨테이너.
`.rds` 형식으로 저장.

```
Seurat Object
│
├── @assays$RNA
│   ├── counts        원본 count matrix
│   ├── data          정규화된 matrix
│   └── scale.data    스케일된 matrix
├── @meta.data        세포 메타데이터 (AnnData .obs 대응)
├── @reductions
│   ├── pca           PCA 임베딩 (AnnData .obsm["X_pca"] 대응)
│   └── umap          UMAP 임베딩 (AnnData .obsm["X_umap"] 대응)
├── @graphs           neighbor graph (AnnData .obsp 대응)
└── @misc             기타 정보 (AnnData .uns 대응)
```

### AnnData ↔ Seurat Object 대응 관계

| AnnData | Seurat Object | 내용 |
|---------|--------------|------|
| `.X` | `@assays$RNA@data` | 정규화된 발현 행렬 |
| `.raw.X` | `@assays$RNA@counts` | 원본 count |
| `.obs` | `@meta.data` | 세포 메타데이터 |
| `.var` | `@assays$RNA@meta.features` | 유전자 메타데이터 |
| `.obsm["X_pca"]` | `@reductions$pca@cell.embeddings` | PCA 좌표 |
| `.obsm["X_umap"]` | `@reductions$umap@cell.embeddings` | UMAP 좌표 |
| `.obsp["connectivities"]` | `@graphs$RNA_nn` | neighbor graph |
| `.uns` | `@misc` | 기타 정보 |

---

## 5. 분석 프레임워크

### Scanpy vs Seurat 기능 비교

| 분석 단계 | Scanpy (Python) | Seurat (R) |
|-----------|----------------|------------|
| QC 지표 계산 | `sc.pp.calculate_qc_metrics()` | `PercentageFeatureSet()` |
| 정규화 | `sc.pp.normalize_total()` | `NormalizeData()` |
| Log 변환 | `sc.pp.log1p()` | `NormalizeData(log=TRUE)` |
| HVG 선정 | `sc.pp.highly_variable_genes()` | `FindVariableFeatures()` |
| 스케일링 | `sc.pp.scale()` | `ScaleData()` |
| PCA | `sc.tl.pca()` | `RunPCA()` |
| Batch correction | `harmonypy` (외부) | `IntegrateLayers()` / `harmony` |
| Neighbor graph | `sc.pp.neighbors()` | `FindNeighbors()` |
| UMAP | `sc.tl.umap()` | `RunUMAP()` |
| 클러스터링 | `sc.tl.leiden()` | `FindClusters()` |
| DEG | `sc.tl.rank_genes_groups()` | `FindMarkers()` |
| Trajectory | `sc.tl.paga()` / scVelo | Monocle3 (별도 패키지) |
| Cell annotation | CellTypist (외부) | `SingleR` (별도 패키지) |

### 프레임워크 선택 기준

```
Python (Scanpy) 선택 권장:
  - 대용량 데이터 (세포 수 > 100,000)
  - 속도/메모리 효율이 중요한 경우
  - 커스텀 파이프라인 자동화 (Snakemake 연동)
  - velocity, spatial 등 최신 분석 도구 접근성

R (Seurat) 선택 권장:
  - 생물학자 중심 팀 (R 친숙도 높음)
  - 레퍼런스/튜토리얼이 많은 분석 방법 필요
  - Bioconductor 생태계 활용 필요
  - 기존 Seurat 코드베이스 유지보수
```

---

## 6. 포맷 변환

### 변환이 필요한 경우

| 상황 | 변환 방향 | 필요 여부 |
|------|----------|----------|
| 같은 생태계 내 분석 | — | ❌ 불필요 |
| Python으로 Seurat 결과 이어받기 | RDS → H5AD | ✅ 필요 |
| R로 Scanpy 결과 이어받기 | H5AD → RDS | ✅ 필요 |
| Monocle3 사용 (trajectory, R 전용) | H5AD → RDS | ✅ 필요 |
| 협업자와 데이터 공유 | 상황에 따라 | 경우에 따라 |

### 변환 방법

**RDS → H5AD (R에서 실행)**

```r
library(Seurat)
library(SeuratDisk)

obj <- readRDS("sample.rds")
SaveH5Seurat(obj, filename = "sample.h5seurat")
Convert("sample.h5seurat", dest = "h5ad")
# → sample.h5ad 생성됨
```

**H5AD → RDS (R에서 실행)**

```r
library(zellkonverter)
adata <- readH5AD("sample.h5ad")          # SingleCellExperiment로 로드
seurat_obj <- as.Seurat(adata)            # Seurat Object로 변환
saveRDS(seurat_obj, "sample_from_scanpy.rds")
```

**CSV + metadata → AnnData (Python에서 실행)**

```python
import scanpy as sc
import pandas as pd

counts = pd.read_csv("Count_data.csv", index_col=0)
metadata = pd.read_csv("Barcode_metadata.csv", index_col=0)

adata = sc.AnnData(X=counts.T)           # genes × cells → cells × genes
adata.obs = metadata.loc[adata.obs_names]
adata.write_h5ad("sample.h5ad")
```

---

## 7. 현재 프로젝트 위치

```
Parse Biosciences 실험
         │
      split-pipe (1차 분석)
         │
      .h5ad 파일 (AnnData)
         │
      Scanpy 파이프라인 ← 현재 구축된 파이프라인
      (Snakemake 자동화)
         │
         ├── QC → 정규화 → HVG → PCA → Harmony → UMAP → Leiden
         │
         └── Jupyter Notebooks
               ├── Cell type annotation (CellTypist)
               ├── DEG 분석
               ├── Trajectory (PAGA + DPT / scVelo)
               └── Figure 생성
```

### 지원 입력 포맷 (`src/modules/io.py`)

| 포맷 | 출처 | 로드 방법 |
|------|------|----------|
| `.h5ad` | split-pipe, alevin, 변환된 Seurat | `sc.read_h5ad()` |
| MTX 디렉토리 | CellRanger, STARsolo | `sc.read_10x_mtx()` |
| `.h5` | CellRanger | `sc.read_10x_h5()` |
| `.loom` | velocyto, kallisto | `sc.read_loom()` |
| `.csv` | Seurat export, 커스텀 | `sc.read_csv()` |

---

## 8. 용어 정리

| 용어 | 정의 |
|------|------|
| **FASTQ** | 시퀀싱 원본 데이터 파일 (reads + quality score) |
| **Count matrix** | 세포 × 유전자 발현량 행렬 (UMI count) |
| **UMI** | Unique Molecular Identifier — PCR 중복 제거용 바코드 |
| **Barcode** | 개별 세포를 구분하는 DNA 서열 |
| **MTX** | Matrix Market Exchange 포맷 — 희소 행렬 저장 표준 |
| **HDF5 / H5** | 계층적 데이터 포맷 — 대용량 데이터 저장에 사용 |
| **H5AD** | AnnData용 HDF5 포맷 (`.h5ad`) |
| **HVG** | Highly Variable Genes — 샘플 간 발현 변이가 큰 유전자 |
| **Doublet** | 하나의 바코드에 두 세포가 잡힌 기술적 오류 |
| **Ambient RNA** | 세포 분리 과정에서 유입된 유리 RNA |
| **Batch effect** | 실험 배치(날짜, 시약 lot 등) 차이로 인한 기술적 변이 |
| **Leiden / Louvain** | 그래프 기반 클러스터링 알고리즘 |
| **PAGA** | Partition-based Graph Abstraction — trajectory 분석법 |
| **RNA velocity** | spliced/unspliced 비율로 분화 방향 추정 |
| **Pseudotime** | 세포의 분화/발달 단계를 수치로 표현한 것 |
