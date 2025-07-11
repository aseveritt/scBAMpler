#!/bin/bash

eval "$(conda shell.bash hook)"

# Create and activate environment
conda create -n macs_Renv pip r-base macs3 zlib -y
conda activate macs_Renv

Rscript helper_scripts/peak_calling/install_R_packages.R

# Optionally, grab species-specific exclusion list
#wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
#gunzip hg38-blacklist.v2.bed.gz
