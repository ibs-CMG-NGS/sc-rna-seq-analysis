"""
src/modules/io.py
AnnData 로드/저장 유틸리티

지원 포맷:
  - h5ad  : AnnData HDF5 (split-pipe, alevin 등)
  - 10x_mtx: 10x Genomics MEX 디렉토리 (CellRanger, STARsolo)
  - 10x_h5 : 10x Genomics HDF5 단일 파일 (CellRanger .h5)
  - loom   : Loom 포맷 (velocyto 출력 등)
  - csv    : 단순 CSV 카운트 매트릭스
  - auto   : 파일 확장자/타입으로 자동 감지
"""

from pathlib import Path
import scanpy as sc
import anndata as ad


_LOADERS = {
    "h5ad": lambda p: sc.read_h5ad(p),
    "10x_mtx": lambda p: sc.read_10x_mtx(p, var_names="gene_symbols", cache=True),
    "10x_h5": lambda p: sc.read_10x_h5(p),
    "loom": lambda p: sc.read_loom(p),
    "csv": lambda p: sc.read_csv(p).T,
}


def load_adata(path: str, fmt: str = "auto") -> ad.AnnData:
    """
    AnnData 로드 (포맷 자동 감지 또는 명시 지정)

    Parameters
    ----------
    path : str
        파일 또는 디렉토리 경로
    fmt : str
        'auto' | 'h5ad' | '10x_mtx' | '10x_h5' | 'loom' | 'csv'

    Returns
    -------
    AnnData
    """
    if fmt == "auto":
        fmt = _detect_format(path)

    if fmt not in _LOADERS:
        raise ValueError(f"Unsupported format: '{fmt}'. Available: {list(_LOADERS)}")

    adata = _LOADERS[fmt](path)
    adata.var_names_make_unique()
    return adata


def save_adata(adata: ad.AnnData, path: str) -> None:
    """AnnData를 h5ad로 저장 (부모 디렉토리 자동 생성)"""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, compression="gzip")


def check_h5ad_has_spliced(data_dir: str) -> bool:
    """
    디렉토리 내 h5ad 파일 중 'spliced' layer를 가진 파일이 있는지 확인
    (RNA velocity를 위한 전처리 여부 감지)
    """
    for h5ad_path in Path(data_dir).rglob("*.h5ad"):
        try:
            adata = sc.read_h5ad(h5ad_path, backed="r")
            if "spliced" in adata.layers:
                return True
        except Exception:
            continue
    return False


def detect_velocity_mode(data_dir: str) -> str:
    """
    RNA velocity 분석 가능 여부 감지

    Returns
    -------
    'scvelo_ready' : loom 파일 또는 spliced layer가 있는 h5ad 발견
    'bam_available': BAM 파일 발견 (velocyto 실행 필요)
    'no_velocity'  : velocity 관련 파일 없음 → PAGA + CytoTRACE 사용
    """
    data_path = Path(data_dir)

    # 1순위: spliced/unspliced가 포함된 loom 파일
    loom_files = list(data_path.rglob("*.loom"))
    if loom_files:
        return "scvelo_ready"

    # 2순위: h5ad에 spliced layer 포함
    if check_h5ad_has_spliced(data_dir):
        return "scvelo_ready"

    # 3순위: BAM 파일 (velocyto 실행 필요)
    bam_files = list(data_path.rglob("*.bam"))
    if bam_files:
        return "bam_available"

    # 없음
    return "no_velocity"


def _detect_format(path: str) -> str:
    p = Path(path)
    if p.suffix == ".h5ad":
        return "h5ad"
    if p.suffix == ".h5":
        return "10x_h5"
    if p.suffix == ".loom":
        return "loom"
    if p.suffix == ".csv":
        return "csv"
    if p.is_dir():
        return "10x_mtx"
    raise ValueError(
        f"Cannot detect format for: {path}\n"
        f"Please specify fmt= explicitly: h5ad | 10x_mtx | 10x_h5 | loom | csv"
    )
