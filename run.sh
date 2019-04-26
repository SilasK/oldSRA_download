snakemake -s SRA_download/Snakefile --config url_table=~/Metagenomics/Data/EODF/EODF_info.tab.txt -d ~/Metagenomics/Data/ --use-conda -j 4 single
