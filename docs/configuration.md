# Configuration Guide

## configs/ 파일 구조

```
configs/
├── pipeline.yaml   # 분석 파라미터 (QC 임계값, 알고리즘 설정 등)
├── samples.yaml    # 입력 데이터 및 샘플 메타데이터
└── markers.yaml    # Cell type 마커 유전자 목록
```

---

## pipeline.yaml 파라미터 설명

### QC 임계값

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `qc.min_genes` | 200 | 세포당 최소 유전자 수 (낮으면 저품질 세포) |
| `qc.max_genes` | 8000 | 세포당 최대 유전자 수 (높으면 doublet 의심) |
| `qc.min_counts` | 500 | 세포당 최소 UMI count |
| `qc.max_pct_mt` | 20 | 미토콘드리아 유전자 비율 최대 (%) |
| `qc.scrublet.expected_doublet_rate` | 0.06 | 예상 doublet 비율 |
| `qc.scrublet.doublet_score_threshold` | 0.25 | doublet 판정 임계값 |

### HVG 선정

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `hvg.flavor` | seurat_v3 | HVG 선정 방법 |
| `hvg.n_top_genes` | 3000 | 선정할 HVG 수 |
| `hvg.batch_key` | sample | batch 보정 기준 칼럼 |

### 클러스터링

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `clustering.resolutions` | [0.3, 0.5, 0.8, 1.0, 1.2] | Leiden 다중 resolution |
| `clustering.default_resolution` | 0.8 | `leiden` 칼럼의 기본 resolution |

---

## samples.yaml 파라미터 설명

### datasets 섹션

| 필드 | 필수 | 설명 |
|------|------|------|
| `path` | ✅ | h5ad 또는 10x 디렉토리 경로 |
| `format` | - | `auto` \| `h5ad` \| `10x_mtx` \| `10x_h5` \| `loom` \| `csv` |
| `species` | - | `human` \| `mouse` \| `mixed` |
| `run_in_pipeline` | - | `true`이면 Snakemake 기본 실행 대상 |

### samples 섹션

샘플 이름은 h5ad의 `obs` 칼럼에 있는 `sample` 값과 일치해야 합니다.

| 필드 | 설명 |
|------|------|
| `condition` | 실험 조건. 확인 전까지 "unknown" 유지 |
| `replicate` | replicate 번호 |

---

## markers.yaml 구조

```yaml
# 분류 체계 (Category)
T_cells:
  # 세포 유형 (Subtype)
  CD4_T:
    markers: ["CD3D", "CD4", "IL7R"]  # 마커 유전자 리스트
    description: "CD4+ T cells"        # 설명 (선택)
```

`annotation.py`의 `score_cell_types()` 및 `plotting.py`의 `marker_dotplot()`에서 이 구조를 사용합니다.
