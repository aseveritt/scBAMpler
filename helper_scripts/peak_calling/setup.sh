#!/bin/bash

eval "$(conda shell.bash hook)"

# Create and activate environment
conda create -n macs_Renv -c conda-forge -c bioconda -y \
    pip r-base macs3 zlib \
    r-optparse r-dplyr r-stringr r-data.table \
    bioconductor-genomicranges \
    bioconductor-genomicfeatures \
    bioconductor-rtracklayer \
    bioconductor-txdb.hsapiens.ucsc.hg38.knowngene

conda activate macs_Renv

# Optionally, grab species-specific exclusion list
#wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
#gunzip hg38-blacklist.v2.bed.gz
