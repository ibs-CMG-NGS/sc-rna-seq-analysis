# Setup Guide

## 요구 사항

- Linux (Ubuntu 20.04+) 또는 WSL2
- Conda (Miniconda3 또는 Anaconda)
- 최소 32GB RAM (74K 세포 처리 기준)
- Python 3.10

---

## 1. Conda 환경 설치

```bash
# 프로젝트 루트에서 실행
conda env create -f environment.yml

# 환경 활성화
conda activate sc-rna-seq-analysis
```

### Jupyter kernel 등록

```bash
conda activate sc-rna-seq-analysis
python -m ipykernel install --user --name sc-rna-seq-analysis --display-name "sc-rna-seq-analysis"
```

---

## 2. 데이터 경로 설정

`configs/samples.yaml`에서 데이터셋 경로를 확인합니다.

현재 등록된 데이터:
```
data/all-sample/DGE_filtered/
  ├── human_ratio95_humanGenesOnly.h5ad   ← 메인 분석 대상
  ├── human_ratio95_cells.h5ad
  ├── human_ratio95_EGFPpos_humanGenesOnly.h5ad
  ├── human_ratio95_EGFPpos_cells.h5ad
  └── trailmaker_raw_mixed.h5ad
```

경로가 다를 경우 `configs/samples.yaml`의 `path` 필드를 수정하세요.

---

## 3. 샘플 조건 메타데이터 확인

`configs/samples.yaml`의 `samples` 섹션에서 각 샘플의 `condition`을 업데이트하세요.

```yaml
samples:
  APC_1:
    condition: "treatment"  # "unknown" → 실제 조건으로 변경
    replicate: 1
```

---

## 4. RNA Velocity 설정 (선택)

BAM 파일 또는 loom 파일이 있는 경우 `configs/samples.yaml`에서:

```yaml
velocity:
  bam_dir: "data/bam/"      # BAM 파일 디렉토리
  loom_dir: "data/loom/"    # loom 파일 디렉토리 (이미 있는 경우)
```

BAM 파일이 있고 scVelo를 사용하려면 추가 설치:
```bash
pip install scvelo velocyto
```

---

## 5. CellTypist 모델 다운로드

annotation notebook 실행 시 자동으로 다운로드되지만, 사전 다운로드 원하면:

```bash
python -c "
import celltypist
from celltypist import models
models.download_models(model='Immune_All_Low.pkl')
"
```

---

## 6. 설치 확인

```bash
conda activate sc-rna-seq-analysis
python -c "
import scanpy as sc
import anndata as ad
import leidenalg
import harmonypy
print('scanpy:', sc.__version__)
print('All packages OK')
"
```
