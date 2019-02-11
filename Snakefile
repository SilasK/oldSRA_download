
import pandas as pd
SRR_list = pd.read_csv(config['url_table'],sep='\t',index_col=0).iloc[:,0]

outdir='.'

rule all:
    input:
        expand("{outdir}/{SRR}_{direction}.fastq.gz",SRR=SRR_list,outdir=outdir,direction=['R1','R2'])


rule download_SRR:
    output:
        "{outdir}/{SRR}_1.fastq.gz",
        "{outdir}/{SRR}_2.fastq.gz"
    params:
        outdir=outdir
    threads:
        4
    conda:
        "envs/download.yaml"
    shell:
        "parallel-fastq-dump --sra-id {wildcards.SRR} --threads {threads} --gzip --split-files --outdir {params.outdir}"

localrules: rename_SRR
rule rename_SRR:
    input:
        "{outdir}/{SRR}_1.fastq.gz",
        "{outdir}/{SRR}_2.fastq.gz"
    output:
        "{outdir}/{SRR}_R1.fastq.gz",
        "{outdir}/{SRR}_R2.fastq.gz"
    threads:
        1
    shell:
        "mv {input[0]} {output[0]}; "
        "mv {input[1]} {output[1]}"
