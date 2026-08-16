#!/bin/bash
# Creates scBAMpler_R_env from scBAMpler_R_env.yaml, the environment for the optional R
# helper scripts (peak calling, and building the extension's H5 input from an ArchR project).
# Not required to run scBAMpler itself.
#
#   bash helper_scripts/setup_scBAMpler_R_env.sh                # everything, including ArchR
#   bash helper_scripts/setup_scBAMpler_R_env.sh --skip-archr   # peak calling only, much smaller

set -eo pipefail

SKIP_ARCHR=0
if [[ "${1:-}" == "--skip-archr" ]]; then SKIP_ARCHR=1; fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

eval "$(conda shell.bash hook)"

conda env create -f "$SCRIPT_DIR/scBAMpler_R_env.yaml"

conda activate scBAMpler_R_env

# Confirm what call_peaks.R needs is importable before the user runs it.
Rscript -e 'suppressPackageStartupMessages({
  library(optparse); library(GenomicRanges); library(GenomicFeatures)
  library(dplyr); library(rtracklayer); library(stringr)
  library(data.table); library(parallel)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
}); cat("Peak calling dependencies loaded successfully.\n")'

if [[ "$SKIP_ARCHR" -eq 0 ]]; then
    # ArchR is not on CRAN, Bioconductor, or conda.
    Rscript -e 'remotes::install_github("GreenleafLab/ArchR", upgrade = "never")'
    Rscript -e 'suppressPackageStartupMessages({
      library(ArchR); library(rhdf5); library(Matrix)
    }); cat("ArchR dependencies loaded successfully.\n")'
fi

echo
echo "Done. Activate the environment before running the helper scripts:"
echo "    conda activate scBAMpler_R_env"

# Optionally, grab species-specific exclusion list
#wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
#gunzip hg38-blacklist.v2.bed.gz
