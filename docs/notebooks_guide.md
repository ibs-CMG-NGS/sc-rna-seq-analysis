# Notebooks Guide

## 실행 전 준비

1. Snakemake 파이프라인 완료 확인:
   ```bash
   ls output/checkpoints/human_genes_only/07_clustered.h5ad
   ```
2. Jupyter kernel 등록:
   ```bash
   python -m ipykernel install --user --name sc-rna-seq-analysis
   ```
3. JupyterLab 실행:
   ```bash
   jupyter lab --notebook-dir=src/notebooks
   ```

---

## Notebook 순서 및 설명

### 01_cluster_review.ipynb

**목적**: Snakemake 클러스터링 결과 검토 및 저품질 클러스터 제거
**입력**: `output/checkpoints/{dataset}/07_clustered.h5ad`
**출력**: `output/checkpoints/{dataset}/07_clustered_clean.h5ad`

**주요 작업**:
- 클러스터별 QC 지표 violin plot 확인
- 높은 mt%, 낮은 gene count 클러스터 식별
- `LOW_QUALITY_CLUSTERS` 리스트에 제거할 클러스터 번호 입력

**설정 변경 위치**:
```python
DATASET = "human_genes_only"          # 분석할 데이터셋
LOW_QUALITY_CLUSTERS = ["12", "25"]   # 제거할 클러스터 (빈 리스트면 제거 안 함)
```

---

### 02_cell_type_annotation.ipynb

**목적**: CellTypist 자동 분류 + 마커 유전자 검증 + 수동 보정
**입력**: `output/checkpoints/{dataset}/07_clustered_clean.h5ad`
**출력**: `output/checkpoints/{dataset}/08_annotated.h5ad`

**주요 작업**:
- CellTypist `Immune_All_Low.pkl` 모델로 자동 분류
- 마커 유전자 dotplot으로 분류 검증
- `cluster_to_celltype` 딕셔너리로 수동 annotation 적용

**설정 변경 위치**:
```python
cluster_to_celltype = {
    "0": "CD4 T",
    "1": "CD8 T",
    # ...
}
```

---

### 03_deg_analysis.ipynb

**목적**: Cell type 간 DEG 분석 (Wilcoxon rank-sum)
**입력**: `output/checkpoints/{dataset}/08_annotated.h5ad`
**출력**: `output/tables/de_results/{dataset}/`, volcano plot

**주요 작업**:
- `sc.tl.rank_genes_groups()` one-vs-rest DEG
- 상위 마커 dotplot
- Volcano plot
- CSV 테이블 export

**조건 간 DEG 활성화**:
`configs/samples.yaml`의 `condition` 채운 후 notebook 하단 주석 해제

---

### 04_trajectory_paga.ipynb

**목적**: PAGA 기반 trajectory 분석 (velocity 데이터 감지 후 자동 분기)
**입력**: `output/checkpoints/{dataset}/08_annotated.h5ad`
**출력**: `output/checkpoints/{dataset}/09_trajectory.h5ad`

**자동 분기 로직**:
- loom 파일 또는 spliced layer 있음 → scVelo
- BAM 파일 있음 → velocyto 실행 안내 후 대기
- 없음 → PAGA + DPT + CytoTRACE2

**설정 변경 위치**:
```python
TRAJECTORY_GENES = ["MKI67", "CD44"]  # pseudotime에 따라 시각화할 유전자
```

---

### 05_figures.ipynb

**목적**: 논문/보고서용 고해상도 figure 생성 (300 DPI, PDF)
**입력**: `output/checkpoints/{dataset}/09_trajectory.h5ad`
**출력**: `output/figures/{dataset}/publication/`

**생성되는 figure**:
- `fig1_umap_overview.pdf/png` — cell type + sample UMAP
- `fig2_marker_dotplot.pdf` — 마커 유전자 dotplot
- `fig3_trajectory.png` — pseudotime/velocity UMAP
- `fig4_cell_proportion.pdf/png` — 샘플별 세포 비율 stacked bar
