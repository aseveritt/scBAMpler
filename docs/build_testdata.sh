#!/usr/bin/env bash

source "$(conda info --base)/etc/profile.d/conda.sh"
set -e

### =====================================================================
### Provenance record for the scBAMpler test data
###
### These are the commands used to produce every file distributed in
### test_data.tar.gz and example_output.tar.gz. You do not need to run any of
### this to use scBAMpler -- the resulting files are on Zenodo. This record
### exists so the inputs can be traced back to their public sources and
### rebuilt from scratch if needed.
###
### All paths assume you are in the repository root.
### =====================================================================


#### -------- Environments -------- ###
# Main pipeline (see README, Installation)

# commented out because im assuming they already exist. 
# conda create -n scBAMpler_env -c conda-forge -c bioconda python=3.10 numpy scipy pandas samtools bedtools sinto "setuptools<81" -y
# conda activate scBAMpler_env
# pip install h5py k-means-constrained    # specific to the cell similarity extension, h5py can be from conda too. 

# # Optional R environment, used for peak calling and for building the H5 input
# bash helper_scripts/setup_scBAMpler_R_env.sh

conda activate scBAMpler_env

#### -------- Retrieve and subsample ENCODE scATAC-seq cell lines -------- ###
# Three cell lines, so the cell-homogeneity extension has multiple populations to mix.
# ENCODE BAMs are coordinate sorted, so sinto's output is too; we only need to index it.
#
#         unfiltered BAM   experiment    cell-call QC (multiplet quality metric UUID)
# HepG2:  ENCFF089OZP      ENCSR398OHC   4381cc94-c5fa-4e65-97cf-3c7b202ce573
# K562:   ENCFF473MZT      ENCSR217VXJ   4458ac03-1e6c-43c1-b206-9adcf1c323b6
# MCF-7:  ENCFF180UYN      ENCSR341TTX   67ea62ae-8dae-41dd-82eb-e05701205ec6

mkdir -p test_data example_output

#keep as for-loop or break into sbatch's on slurm
for LINE in HEPG2 K562 MCF7; do
    case $LINE in
        HEPG2) ACC=ENCFF089OZP ; QC=4381cc94-c5fa-4e65-97cf-3c7b202ce573 ;;
        K562)  ACC=ENCFF473MZT ; QC=4458ac03-1e6c-43c1-b206-9adcf1c323b6 ;;
        MCF7)  ACC=ENCFF180UYN ; QC=67ea62ae-8dae-41dd-82eb-e05701205ec6 ;;
    esac

    #download bams
    wget https://www.encodeproject.org/files/${ACC}/@@download/${ACC}.bam -O test_data/${LINE}_unfiltered.bam

    # get ENCODE's per-barcode cell calls, attached as QC to each experiment's fragments file.
    # Barcodes the pipeline called as cells have IsMultiplet == False; the rest are "Indeterminate" (background).
    wget https://www.encodeproject.org/sc-atac-multiplet-quality-metrics/${QC}/@@download/barcodes_status/multiplet_barcodes_status.tsv \
        -O test_data/${LINE}_barcodes_status.tsv

    # keep ENCODE-called cells, then a random 2,000 of them
    awk -F'\t' '$2=="False" {print $1}' test_data/${LINE}_barcodes_status.tsv > test_data/${LINE}_barcodes_all.txt
    shuf -n 2000 test_data/${LINE}_barcodes_all.txt | awk -v OFS="\t" '{print $1,"'"${LINE}"'_subset"}' > test_data/${LINE}_barcodes_sub.tsv

    # the label in column 2 determines the output filename, i.e. <LINE>_subset.bam
    sinto filterbarcodes \
        --bam test_data/${LINE}_unfiltered.bam \
        --cells test_data/${LINE}_barcodes_sub.tsv \
        --outdir test_data/ \
        --nproc 8

    samtools quickcheck -v test_data/${LINE}_subset.bam
    samtools index test_data/${LINE}_subset.bam
done

#rm test_data/*_barcodes_sub.tsv test_data/_barcodes_status.tsv test_data/*_unfiltered.bam*

#### -------- Exclusion list -------- ###
# The distributed peak files were called with the ENCODE hg38 blacklist.
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz -O test_data/hg38-blacklist.v2.bed.gz
gunzip test_data/hg38-blacklist.v2.bed.gz


#### -------- Call Peaks -------- ###
conda activate scBAMpler_R_env

for LINE in HEPG2 K562 MCF7; do
    Rscript helper_scripts/peak_calling/call_peaks.R \
        --bam_file test_data/${LINE}_subset.bam \
        --outdir test_data/ \
        --peak_length 500 \
        --exclusion_file test_data/hg38-blacklist.v2.bed \
        --cores 8
done

# Union peak set across all three lines. This defines the row order of the peak matrix in peakmat_input.h5.
Rscript helper_scripts/peak_calling/call_peaks.R \
    --outdir test_data/ \
    --union_files "HEPG2_subset_standardized_500bp.bed,K562_subset_standardized_500bp.bed,MCF7_subset_standardized_500bp.bed" \
    --union_outfile union_standardized_500bp.bed \
    --cores 8


#### -------- Build the H5 input for the extension -------- ###
# Requires ArchR, installed by setup_scBAMpler_R_env.sh above.
# Input paths and sample names are constants at the top of the script.
Rscript helper_scripts/H5_from_ArchR/MakeH5.R


#### -------- example_output.tar.gz contents -------- ###
conda activate scBAMpler_env

scBAMpler create-dictionary \
    --bam_file test_data/HEPG2_subset.bam \
    --peak_file test_data/HEPG2_subset_standardized_500bp.bed \
    --output_file example_output/HEPG2_subset.pickle \
    --verbose

scBAMpler sampler \
    --input_pickle example_output/HEPG2_subset.pickle \
    --input_bam test_data/HEPG2_subset.bam \
    --output_prefix example_output/HEPG2_subset_c500_s12 \
    --downsample_by cells \
    --downsample_to 500 \
    --seed 12 \
    --nproc 10 \
    --output_fragment \
    --verbose

scBAMpler make-pseudobulks \
    --input test_data/peakmat_input.h5 \
    --output example_output/medoids_s50.pickle \
    --dimred umap \
    --label-col CellLine \
    --cluster-size 50 \
    --nproc 5

scBAMpler mix-pseudobulks \
    --input example_output/medoids_s50.pickle \
    --output example_output/combos_all.csv \
    --groups all \
    --n-combos 2000 \
    --cluster-size 50 \
    --ft-sizes 200 300 400 500 600 700 800

scBAMpler mix-pseudobulks \
    --input example_output/medoids_s50.pickle \
    --output example_output/combos_k562_hepg2.csv \
    --groups K562 HEPG2 \
    --n-combos 500 \
    --cluster-size 50 \
    --ft-sizes 200 250 300 350 400

cat example_output/combos_all.csv <(tail -n +2 example_output/combos_k562_hepg2.csv) > example_output/combos_combined.csv


#run notebook. 

scBAMpler extract-populations \
    --barcode-dir example_output/selected/barcodes \
    --output example_output/populations/ \
    --bam HEPG2=test_data/HEPG2_subset.bam \
    --bam K562=test_data/K562_subset.bam \
    --bam MCF7=test_data/MCF7_subset.bam \
    --nproc 8

bash example_output/populations/scripts/combo_0.sh



