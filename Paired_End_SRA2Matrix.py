#######################################
### For Bulk RNA seq Snakemake file ####
### From SRA to count matrix ####
#### For paied end file ####
##### Boyuan Wang 20230405 ####
#########################################
import os
data_dir="SRA"
SAMPLES=os.listdir(data_dir)
i=0
while i<len(SAMPLES):
        SAMPLES[i]=SAMPLES[i].split('.')[0]
        i+=1
samples=list(set(SAMPLES))
INDEX_HISAT2_mm39="/DATA/public/reference/hisat2/mm39/genome"
GTF_mm39="/DATA/public/reference/GENCODE/mm39/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf"
rule all:
    input:
        "qc/multiqc_cutadapter.html",
        "counts.txt",
        "counts.txt.summary"
rule fastqdump:
    input:
         "SRA/{sample}.sra"
    output:
         "fastq/{sample}_1.fastq",
         "fastq/{sample}_2.fastq",
    log:
        "logs/{sample}.SRA2fastq.log"
    shell:
        "/opt/mambaforge/envs/rna/bin/fastq-dump  --split-files -O fastq/ {input} >{log} 2>&1"
rule pigz:
    input:
        "fastq/{sample}_1.fastq",
        "fastq/{sample}_2.fastq",
    output:
        "fastq/{sample}_1.fastq.gz",
        "fastq/{sample}_2.fastq.gz",
    log:
        "logs/{sample}.fastq2fqgz.log"
    shell:
        "pigz -p 20 {input} >{log} 2>&1"
rule trim_galore:
    input:
        ["fastq/{sample}_1.fastq.gz", "fastq/{sample}_2.fastq.gz"],
    output:
        "trimmed/{sample}_1_val_1.fq.gz",
        "trimmed/{sample}_1.fastq.gz_trimming_report.txt",
        "trimmed/{sample}_2_val_2.fq.gz",
        "trimmed/{sample}_2.fastq.gz_trimming_report.txt",
    log:
        "logs/trim_galore.{sample}.log",
    shell:
        "trim_galore  -j 4 -q 25 --phred33 --length 36 -paired -o trimmed/ {input} > {log} 2>&1"
rule multiqc_cutadapter:
    input:
        expand("trimmed/{sample}_{read}.fastq.gz_trimming_report.txt", sample=samples,read=["1","2"])
    output:
        "qc/multiqc_cutadapter.html"
    log:
        "logs/multiqc_cutadapter.log"
    shell:
        "multiqc -f -n multiqc_cutadapter.html -o qc trimmed > {log} 2>&1"
rule hisat2:
    input:
        "trimmed/{sample}_1_val_1.fq.gz",
        "trimmed/{sample}_2_val_2.fq.gz"
    output:
        "sam/{sample}.sam"
    log:
        "logs/hisat2.{sample}.log"
    shell:
        "hisat2 -x {INDEX_HISAT2_mm39} -p 10 -1 {input[0]} -2 {input[1]} -S {output} > {log} 2>&1"
rule samtools_sort:
    input:
        "sam/{sample}.sam"
    output:
        "bam_sorted/{sample}_sorted.bam"
    log:
        "logs/samtools_sort.{sample}.log"
    shell:
        "samtools sort -n -@ 15 {input} -o {output} >{log} 2>&1"
rule featureCounts:
    input:
        expand("bam_sorted/{sample}_sorted.bam",sample=samples)
    output:
        "counts.txt",
        "counts.txt.summary"
    log:
        "logs/featureCounts.log"
    shell:
        "featureCounts -p -B -t exon -g gene_name -a {GTF_mm39} -T 5 -o {output[0]} {input} > {log} 2>&1"
