
import pandas as pd
SRR_list = "SRR8291352 SRR5927703 SRR8291346 SRR8291368 SRR8291348 SRR8291350 SRR8291347 SRR8291351 SRR8291343 SRR8291367 SRR5580763 SRR5927721".split()#"SRR6216944 SRR7761329 SRR8291345 SRR5927712".split()
#
rule all:
    input:
        expand("{SRR}/binning/DASTool/checkm/taxonomy.tsv",SRR=SRR_list)

path='/home/kiesers/scratch/NewSamples'
sample_dir= path+'/samples'
sample_table= 'samples.tsv'

def raw_files(wildcards):
    return expand("{outdir}/{SRR}_{direction}.fastq.gz",SRR=wildcards.SRR,outdir=sample_dir,direction=['R1','R2'])
localrules: add_to_sampletable
rule add_to_sampletable:
    input: raw_files
    output: touch("added_{SRR}_to_sampletable")

    run:
        ST= pd.read_csv(sample_table,sep='\t',index_col=0)
        ST.loc[wildcards.SRR,['Reads_raw_R1' ,   'Reads_raw_R2' ,   'BinGroup']] = (*input,wildcards.SRR)
        ST.to_csv(sample_table,sep='\t')

localrules: atlas_binning, atlas_qc,atlas_assembly
rule atlas_qc:
    input:
        raw_files,
        rules.add_to_sampletable.output,
    output:
        "{SRR}/sequence_quality_control/finished_QC"
    threads:
        1
    log:
        "logs/Pipline/QC/{SRR}.log"
    shell:
        """

        atlas run None "{output}"  \
         --profile EL7 --rerun-incomplete \
        -j 3  &> {log}

        """


rule atlas_assembly:
    input:
        qc=rules.atlas_qc.output,
        raw= raw_files
    output:
        "{SRR}/finished_assembly"
    threads:
        1
    log:
        "logs/Pipline/assembly/{SRR}.log"
    shell:
        """

        atlas run None "{output}"   \
         --profile EL7  \
         --restart-times=3 --rerun-incomplete \
        -j 3 &> {log}

        """



rule atlas_binning:
    input:
        qc=rules.atlas_qc.output,
        assembly= rules.atlas_assembly.output
    output:
        directory("{SRR}/binning/DASTool/bins"),
        "{SRR}/binning/DASTool/checkm/taxonomy.tsv"
    threads:
        1
    log:
        "logs/Pipline/binning/{SRR}.log"
    shell:
        """



        atlas run None "{output}"   \
         --profile EL7 --rerun-incomplete \
        -j 5 --restart-times=2 &> {log}

        """


rule download_SRR_paired:
    output:
        "{outdir}/{SRR}_1.fastq.gz",
        "{outdir}/{SRR}_2.fastq.gz"

    #wildcard_constraints:
    #    SRR="SRR[A-Z0-9]+"
    threads:
        4
    conda:
        "envs/download.yaml"
    shell:
        "parallel-fastq-dump --sra-id {wildcards.SRR} --threads {threads} --gzip --split-files --outdir {wildcards.outdir}"

localrules: rename_SRR
rule rename_SRR:
    input:
        "{outdir}/{SRR}_1.fastq.gz",
        "{outdir}/{SRR}_2.fastq.gz"
    output:
        temp("{outdir}/{SRR}_R1.fastq.gz"),
        temp("{outdir}/{SRR}_R2.fastq.gz")
    threads:
        1
    shell:
        "mv {input[0]} {output[0]}; "
        "mv {input[1]} {output[1]}"
