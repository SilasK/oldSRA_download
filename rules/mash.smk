import os

input_dir="MAGs_DasTool2"
Merged_sketch=input_dir+'.msh'
extension="fasta" #"fna.gz"
input_path=os.path.join(input_dir,"{genome}."+extension)

Genomes = glob_wildcards(input_path).genome


rule all:
    input:
        Merged_sketch

rule sketch:
    input:
        input_path
    output:
        "sketches/{genome}.msh"
    params:
        prefix= lambda wc, output: os.path.splitext(output[0])[0]
    threads:
        1
    shell:
        "mash sketch -p {threads} -o {params.prefix} {input}"$




rule paste:
    input:
        expand(rules.sketch.output[0],genome=Genomes)
    output:
        Merged_sketch
    params:
        prefix= lambda wc, output: os.path.splitext(output[0])[0]
    shell:
        " mash paste  {params.prefix} {input}"
