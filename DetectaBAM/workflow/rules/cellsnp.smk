from pathlib import Path

rule cellsnp:
    input:  
        bam="/scratch/nikoe00/bam_dedup/{line}_pseudobulk_unique.dedup.bam",
        barcodes=lambda wc: LINES.at[wc.line, "barcodes"],
        mutations="resources/SomaticMutations/{line}/{line}_SomaticMutations_all.vcf"
    output:
        "results/cellsnp/{line}/cellSNP.base.vcf.gz",
    conda:
        "../envs/cellsnp.yaml",
    log:  
        "logs/cellsnp/{line}.log",
    params:
        out_dir="results/cellsnp/{line}",
    threads: 4
    shell:
        "cellsnp-lite "
        "-s {input.bam} "
        "-b {input.barcodes} "
        "-R {input.mutations} "
        "-O {params.out_dir} "
        "--chrom chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX "
        "--cellTAG CB --UMItag UB "
        "--minCOUNT 1 --minMAF 0 "
        "-p {threads} --gzip "
        "&> {log}"

# rule callability_summary:
#     input:  
#         expand("results/callability/{line}_depth.tsv", line=LINES.index)
#     output: 
#         "results/callability/summary.tsv"
#     conda: 
#         "../envs/samtools.yaml"
#     params: 
#         thr=config["thresholds"]
#     script:
#         "../scripts/summarize.py"