
import pandas as pd
SRR_list = pd.read_csv(config['url_table'],sep='\t',index_col=0).index


if 'outdir' in config:
    outdir = config['outdir']
else:
    outdir=   config['url_table'].replace("_info.tab.txt",'')



rule paired:
    input:
        expand("{outdir}/{SRR}_{direction}.fastq.gz",SRR=SRR_list,outdir=outdir,direction=['R1','R2']),
        #expand("{outdir}/{SRR}.msh",SRR=SRR_list,outdir=outdir)

rule single:
    input:
        expand("{outdir}/{SRR}.fastq.gz",SRR=SRR_list,outdir=outdir),

rule download_SRR_single:
    output:
        "{outdir}/{SRR}.fastq.gz",
    wildcard_constraints:
        SRR="E[A-Z0-9]+"
    params:
        outdir=outdir
    threads:
        4
    conda:
        "envs/download.yaml"
    shell:
        "parallel-fastq-dump --sra-id {wildcards.SRR} --threads {threads} --gzip --outdir {params.outdir}"



rule download_SRR_paired:
    output:
        "{outdir}/{SRR}_1.fastq.gz",
        "{outdir}/{SRR}_2.fastq.gz"
    params:
        outdir=outdir
    wildcard_constraints:
        SRR="E[A-Z0-9]+"
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


rule sketch_reads:
    input:
        "{folder}/{sample}_R1.fastq.gz",
        "{folder}/{sample}_R2.fastq.gz"
    output:
        "{folder}/{sample}.msh"
    params:
        prefix= lambda wc, output: os.path.splitext(output[0])[0]
    conda:
        "envs/mash.yaml"
    threads:
        4
    shell:
        "mash sketch -p {threads} -o {params.prefix} {input}"
