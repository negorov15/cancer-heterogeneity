def pool_bam(wc):  return f"/scratch/nikoe00/bam_data/{LINES.loc[wc.line, 'pool']}/outs/final.bam"

rule subset:
    input:  
        bam=pool_bam, 
        bai=lambda wc: pool_bam(wc)+".bai",
        barcodes=lambda wc: LINES.at[wc.line, "barcodes"]
    output: 
        "/scratch/nikoe00/bam_subsets/{line}_pseudobulk.bam"
    conda: 
        "../envs/subset.yaml"
    log:   
        "logs/subset/{line}.log"
    threads:
        2
    shell:
        "(subset-bam --bam {input.bam} --cell-barcodes {input.barcodes} "
        "--out-bam {output} --cores {threads} && "
        "test $(samtools view -c {output}) -gt 0) > {log} 2>&1"

rule index:
    input:  
        "{prefix}.bam"
    output: 
        "{prefix}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    threads:
        2
    shell:  
        "samtools index -@ {threads} {input}"

rule dedup:
    input:  
        bam="/scratch/nikoe00/bam_subsets/{line}_pseudobulk.bam", 
        bai="/scratch/nikoe00/bam_subsets/{line}_pseudobulk.bam.bai"
    output: 
        "/scratch/nikoe00/bam_dedup/{line}_pseudobulk_unique.dedup.bam"
    conda:
        "../envs/umi_tools.yaml"
    log:   
        "logs/dedup/{line}.log"
    threads:
        1
    shell:
        "umi_tools dedup --stdin={input.bam} --stdout={output} --log={log} --method=unique "
        "--extract-umi-method=tag --umi-tag=UB --cell-tag=CB "
        "> {log} 2>&1"

rule bamcoverage:
    input:  
        bam="/scratch/nikoe00/bam_dedup/{line}_pseudobulk_unique.dedup.bam", 
        bai="/scratch/nikoe00/bam_dedup/{line}_pseudobulk_unique.dedup.bam.bai"
    output:
        "results/bw/{line}.bw"
    conda:
        "../envs/deeptools.yaml"
    log: 
        "logs/bamcoverage/{line}.log"
    threads:
        4
    shell:  
        "bamCoverage --bam {input.bam} -o {output} --binSize 1 --normalizeUsing CPM -p {threads} > {log} 2>&1"

rule compute_matrix:
    input:  
        bw="results/bw/{line}.bw", 
        pos=lambda wc: LINES.at[wc.line, "mutated_genes"]
    output: 
        "results/matrix/{line}_tes_matrix.gz"
    conda: 
        "../envs/deeptools.yaml"
    log:
        "logs/compute_matrix/{line}.log"
    threads:
        4
    shell:
        "computeMatrix --referencePoint TES -R {input.bed} "
        "-S {input.bw} -o {output} -b 3000 -a 500 -p {threads} > {log} 2>&1"