snakemake -s SRA_download/Snakefile --config url_table=/home/kiesers/scratch/miBCgenomes_info.tab.txt -d ~/Metagenomics/Data/ --use-conda -j 4 paired $@
