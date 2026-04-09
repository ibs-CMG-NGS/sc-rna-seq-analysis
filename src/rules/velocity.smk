# src/rules/velocity.smk
# velocyto 실행 → per-sample loom 생성 → scVelo 준비
#
# 전제 조건:
#   configs/pipeline.yaml 의 velocity.gtf 에 GTF 파일 경로 설정
#   configs/samples.yaml 의 velocity.bam_dir 에 BAM 상위 디렉토리 경로 설정
#
# BAM 파일 구조 (Parse Biosciences Split-pipe 출력):
#   {bam_dir}/output_{sample}/barcode_headAligned_anno.sorted.bam
#   예) /data/run1/output_APC_1/barcode_headAligned_anno.sorted.bam
#
# 실행 방법:
#   snakemake --snakefile src/Snakefile --cores 8 \
#     output/velocity/loom/human_genes_only/merged.loom

from pathlib import Path

_vel_cfg = config.get("velocity", {})
_bam_dir  = _vel_cfg.get("bam_dir", "")
_loom_dir = _vel_cfg.get("loom_dir", "")
_gtf      = _vel_cfg.get("gtf", "")
_bam_name = _vel_cfg.get("bam_filename", "barcode_headAligned_anno.sorted.bam")

# velocity 룰은 bam_dir + gtf 가 설정된 경우에만 활성화
if _bam_dir and _gtf:

    rule bam_index:
        """BAI 인덱스 생성 (BAM은 이미 sorted, velocyto 필수 요건)"""
        input:
            bam = str(Path(_bam_dir) / "output_{sample}" / _bam_name),
        output:
            bai = str(Path(_bam_dir) / "output_{sample}" / (_bam_name + ".bai")),
        shell:
            "samtools index {input.bam}"

    rule extract_barcodes:
        """h5ad에서 샘플별 바코드 TSV 추출 (velocyto -b 옵션 입력)"""
        input:
            h5ad = "output/checkpoints/{dataset}/08_annotated.h5ad",
        output:
            barcodes = "output/velocity/barcodes/{dataset}/{sample}.tsv",
        run:
            import scanpy as sc
            from pathlib import Path
            adata = sc.read_h5ad(input.h5ad)
            if "sample" in adata.obs.columns:
                mask = adata.obs["sample"] == wildcards.sample
                barcodes = adata.obs_names[mask].tolist()
            else:
                barcodes = adata.obs_names.tolist()
            Path(output.barcodes).parent.mkdir(parents=True, exist_ok=True)
            with open(output.barcodes, "w") as f:
                f.write("\n".join(barcodes) + "\n")
            print(f"Extracted {len(barcodes)} barcodes for {wildcards.sample}")

    rule velocyto_run:
        """velocyto run: BAM → loom (spliced / unspliced / ambiguous counts)"""
        input:
            bam      = str(Path(_bam_dir) / "output_{sample}" / _bam_name),
            bai      = str(Path(_bam_dir) / "output_{sample}" / (_bam_name + ".bai")),
            barcodes = "output/velocity/barcodes/{dataset}/{sample}.tsv",
            gtf      = _gtf,
        output:
            loom = "output/velocity/loom/{dataset}/{sample}.loom",
        params:
            tmp_dir  = "output/velocity/loom/{dataset}/{sample}_tmp",
            rmsk_flag = lambda wc: (
                f"-m {config['velocity']['repeat_mask_gtf']}"
                if config.get("velocity", {}).get("repeat_mask_gtf")
                else ""
            ),
        shell:
            """
            mkdir -p {params.tmp_dir}
            velocyto run \
                -b {input.barcodes} \
                -o {params.tmp_dir}/ \
                {params.rmsk_flag} \
                {input.bam} \
                {input.gtf}
            mv {params.tmp_dir}/*.loom {output.loom}
            rm -rf {params.tmp_dir}
            """

    rule merge_loom:
        """샘플별 loom 파일을 하나로 병합 → scVelo 입력"""
        input:
            looms = expand(
                "output/velocity/loom/{{dataset}}/{sample}.loom",
                sample=list(config.get("samples", {}).keys()),
            ),
        output:
            merged = "output/velocity/loom/{dataset}/merged.loom",
        run:
            import loompy
            from pathlib import Path
            loompy.combine(list(input.looms), output.merged, key="Accession")
            print(f"Merged {len(input.looms)} loom files → {output.merged}")
