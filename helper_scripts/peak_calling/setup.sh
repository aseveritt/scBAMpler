#!/bin/bash

eval "$(conda shell.bash hook)"

# Create and activate environment
conda env create -f helper_scripts/peak_calling/macs_env.yml
conda activate macs_Renv

Rscript install_R_packages.R

# Optionally, grab species-specific exclusion list
#wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
#gunzip hg38-blacklist.v2.bed.gz
