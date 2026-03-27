# Roadmap

## 현재 구현됨 (v1.0)

### Snakemake Pipeline
- [x] QC 지표 계산 (`compute_qc`)
- [x] 세포 필터링 + Scrublet doublet 제거 (`filter_cells`)
- [x] 정규화 + log1p (`normalize`)
- [x] HVG 선정 + 스케일링 (`select_hvg`)
- [x] PCA (`pca`)
- [x] Harmony batch correction + neighbors + UMAP (`batch_correct_umap`)
- [x] Leiden 클러스터링 (다중 resolution) (`cluster`)

### Jupyter Notebooks
- [x] 클러스터 QC 검토 및 저품질 제거 (`01_cluster_review`)
- [x] Cell type annotation — CellTypist + 수동 (`02_cell_type_annotation`)
- [x] DEG 분석 — Wilcoxon, volcano plot (`03_deg_analysis`)
- [x] Trajectory — PAGA + DPT + CytoTRACE / scVelo 분기 (`04_trajectory_paga`)
- [x] 출판용 figure 생성 (`05_figures`)

### 지원 입력 포맷
- [x] h5ad (split-pipe, alevin 등)
- [x] 10x MTX (CellRanger, STARsolo)
- [x] 10x H5 (CellRanger .h5)
- [x] loom (velocyto)

---

## 계획된 기능 (v1.x)

### 전처리 개선
- [ ] **Ambient RNA 제거 (DecontX)** — snRNA-seq 데이터 시 중요
  - 우선순위: 높음
  - 추가 위치: `src/rules/decontx.smk`, `filter_cells` rule 앞
  - 패키지: `decoupler` 또는 R `celda` 패키지 래퍼

- [ ] **다중 batch correction 옵션** — scVI, BBKNN
  - 현재: Harmony만 지원
  - 추가 위치: `dim_reduction.smk`에 옵션 파라미터

### 추가 분석
- [ ] **Cell-cell communication** (LIANA / CellChat)
  - notebook: `06_cell_communication.ipynb`
  - 패키지: `pip install liana`

- [ ] **Cell composition 통계** (scCODA)
  - 조건 간 세포 유형 비율 변화의 통계적 검정
  - notebook: `07_cell_composition.ipynb`
  - 패키지: `pip install sccoda`

- [ ] **Regulon 분석** (pySCENIC)
  - TF 조절 네트워크
  - notebook: `08_regulon_scenic.ipynb`
  - 패키지: `pip install pyscenic`

- [ ] **Pathway activity** (decoupler)
  - 기본 환경에 `decoupler` 포함됨 → notebook만 추가
  - notebook: `09_pathway_activity.ipynb`

### 인프라
- [ ] **Snakemake profile** — SLURM/SGE HPC cluster 실행 지원
- [ ] **HTML 리포트 자동 생성** — `snakemake --report`
- [ ] **pytest 기반 테스트** — 소규모 샘플 데이터로 파이프라인 검증

---

## 검토 중

- Copy number variation (InferCNV, CopyKAT) — 암 샘플일 경우
- Multi-omics 통합 (muon — CITE-seq, ATAC-seq)
- Cross-dataset integration (scANVI, Seurat v5)
- Label transfer from reference atlas

---

## 변경 이력

| 버전 | 날짜 | 내용 |
|------|------|------|
| v1.0 | 2026-03-27 | 초기 파이프라인 구축 |
