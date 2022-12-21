import os
data_dir="SRA"
SAMPLES=os.listdir(data_dir)
i=0
while i<len(SAMPLES):
        SAMPLES[i]=SAMPLES[i].split('.')[0]
        i+=1
samples=list(set(SAMPLES))
rule all:
    input:    
        expand("fastq/{sample}_S1_L001_I1_001.fastq.gz", sample=samples),
        expand("fastq/{sample}_S1_L001_R1_001.fastq.gz", sample=samples),
        expand("fastq/{sample}_S1_L001_R2_001.fastq.gz", sample=samples)
rule fastqdump:
    input:
         "SRA/{sample}.sra"
    output:
         "fastq/{sample}_1.fastq",
         "fastq/{sample}_2.fastq",
         "fastq/{sample}_3.fastq"
    log:
        "logs/{sample}.SRA2fastq.log"
    shell:
        "/home/wangby/.conda/envs/wby/bin/fastq-dump --split-files -O fastq/ {input} >{log} 2>&1"
rule pigz:
    input:
        "fastq/{sample}_1.fastq",
        "fastq/{sample}_2.fastq",
        "fastq/{sample}_3.fastq"
    output:
        "fastq/{sample}_1.fastq.gz",
        "fastq/{sample}_2.fastq.gz",
        "fastq/{sample}_3.fastq.gz"
    log:
        "logs/{sample}.fastq2fqgz.log"
    shell:
        "pigz -p 20 {input} >{log} 2>&1"
rule changenames:
    input:
        "fastq/{sample}_1.fastq.gz",
        "fastq/{sample}_2.fastq.gz",
        "fastq/{sample}_3.fastq.gz"
    output:
        "fastq/{sample}_S1_L001_I1_001.fastq.gz",
        "fastq/{sample}_S1_L001_R1_001.fastq.gz",
        "fastq/{sample}_S1_L001_R2_001.fastq.gz"
    shell:
        """
        mv {input[0]} {output[0]}
        mv {input[1]} {output[1]}
        mv {input[2]} {output[2]}
        """
