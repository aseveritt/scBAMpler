
<p align ="center">
<img src='https://github.com/aseveritt/scBAMpler/blob/main/docs/scBAMpler.png' style="max-width: 100%; height: auto;">
</p>

scBAMpler was developed to alter one aspect of a scATAC-seq dataset at a time: read count, cell count, and fraction of reads in peaks (FRiP) while preserving the original cell attributes. An extension was further developed which, using multiple cell types, will alter the cell-to-cell homogeneity of a population. For details, please see:

[Comparative evaluation of genomic footprinting algorithms for predicting transcription factor binding sites in single-cell data.](TBD)

## Installation

Please clone the repository:

    $ git clone git+https://github.com/aseveritt/scBAMpler.git

Then, create an environment with required dependencies. Installation and information about miniforge can be found [here](https://github.com/conda-forge/miniforge)

    $ conda create -n scBAMpler_env python=3.10 numpy scipy pandas samtools bedtools sinto -y
    $ conda activate scBAMpler_env
    # cd scBAMpler/
    $ pip install .

---------------

## Download Test Data
```
$ cd scBAMpler
$ scp -R ### test_data/
#roughly XXGb in size

# optionally, if you would like to see the example output directory without running the tutorial:
$ scp -R ## ./example_output
#roughly XXGb in size
```

---------------

## Data Quality Usage

### 1. Call Peak Locations

First, prepare your peak file. This can be done using any method you prefer—the only strict requirement is that the file be in **BED6 format**.
If you would like to use the peak standardization code from our manuscript, we provide it here:

```
$ bash helper_scripts/peak_calling/setup.sh
# ~15 min if fresh-R environment

$ Rscript helper_scripts/peak_calling/call_peaks.R \
    --bam_file test_data/HEPG2_subset.bam \
    --outdir test_data/ \
    --peak_length 500 \
    --cores 8

# ~6 min on subset (2.8Gb), ~35 min on full set (25Gb)
```

#### Input Parameters
* `--bam_file`  
    - Path to the coordinate-sorted input BAM file.
* `--outdir`  
    - Directory where output file will be saved.
* `--peak_length`  
    - Length to which all peaks will be standardized.
* `--txdb`  
    - TxDb package used to get chromosome lengths.  
      Default: `TxDb.Hsapiens.UCSC.hg38.knownGene`
* `--cores` *(optional)*  
    - Number of cores to use with `mclapply`.
* `--exclusion_file` *(optional)*  
    - BED file listing regions to exclude from peak calling.
* `--summit_file` *(optional)*  
    - Use this if a MACS3 file already exists to run only the standardization step.

#### Output
* `<outdir>/*_standardized_<peak_length>bp.bed`  
    - Standardized peaks in BED6 format.





### 2. Build Cell Type Input Dictionaries
Next, build a dictionary for each cell type you want to downsample.  
We assume the BAM file contains a cell barcode tag in the form `CB:Z:*`.

```
$ scBAMpler create-dictionary \
    --bam test_data/HEPG2_subset.bam \
    --peak_file test_data/HEPG2_subset_standardized_500bp.bed \
    --output_file example_output/HEPG2_subset.pickle \
    --verbose

# ~9 min on subset (2.8Gb), ~50 min on full set (25Gb)
```
   
#### Input parameters  
* `--bam`
    - Path to the coordinate-sorted input BAM file.
* `--peak_file`
    - Path to the peak file in BED6 format.
* `--out_file`
    - Path where the final dictionary will be saved (as a `.pickle` file).
* `--verbose`
    - Prints additional progress messages

#### Output    
* `<output_file>`  
    - A Python pickle file containing a dictionary of all cell barcodes, their mapping to peak and non-peak reads, and the necessary numeric encoders. (e.g. HEPG2_subset.pickle)
* `<outfile>`.summary.txt
    - A plain-text file with summary statistics about the cell type.  (e.g. HEPG2_subset.summary.txt)
* `<outfile>`.reads_in_peaks.bed.gz
    - A gzipped BED-like file with two columns: cell barcode and associated peak-read QNAMEs.  
      *(Optionally deleted using the `--delete_intersect` flag.)*
    



### 3. Strategically Downsample BAM
Specify which feature to downsample and to what extent.  
Maximum values are roughly outlined in the `.summary.txt` file generated in the previous step.  
For FRiP, these limits are harder to estimate, but the program will give warning if the requested FRiP is considered too extreme.

```
$ scBAMpler sampler \
    --input_pickle example_output/HEPG2_subset.pickle \
    --input_bam test_data/HEPG2_subset.bam \
    --output_prefix example_output/HEPG2_subset_c500_s12 \
    --downsample_by cells \
    --downsample_to 500 \
    --seed 12 \
    --nproc 10 \
    --output_fragment \
    --verbose

$ scBAMpler sampler \
    --input_pickle example_output/HEPG2_subset.pickle \
    --input_bam test_data/HEPG2_subset.bam \
    --output_prefix example_output/HEPG2_subset_r1e6_s45 \
    --downsample_by reads \
    --downsample_to 1000000 \
    --seed 45 \
    --nproc 10 \
    --output_fragment \
    --verbose

$ scBAMpler sampler \
    --input_pickle example_output/HEPG2_subset.pickle \
    --input_bam test_data/HEPG2_subset.bam \
    --output_prefix example_output/HEPG2_subset_f0.2_s33 \
    --downsample_by frip \
    --downsample_to 0.2 \
    --seed 33 \
    --nproc 10 \
    --output_fragment \
    --verbose

# ~4 min on subset (2.8Gb), ~XX min on full set (25Gb)
```

#### Input Parameters
* `--input_pickle`  
    - Path to the pickle file generated in Step 2.
* `--input_bam`  
    - Path to the coordinate-sorted input BAM file. Reads will be directly extracted from this file.
* `--output_prefix`  
    - Prefix for all output files.
* `--downsample_by`  
    - Type of downsampling to perform.  
      Choices: `"cells"`, `"reads"`, or `"frip"`.
* `--downsample_to`  
    - Target value for the downsampling operation.  
      For `frip` value should be between 0 and 1 (e.g., 0.2) and this represents the pseudobulk-level FRiP (not average per-cell)
* `--seed`  
    - Random seed for reproducibility.
* `--nproc`  
    - Number of processors to use.
* `--output_fragment`  
    - If set, will also output a `fragment.tsv.bgz` file in addition to the BAM file.
* `--verbose`
    - Prints additional progress messages

#### Output
* `<output_prefix>.bam`  
    - Downsampled, coordinate-sorted BAM file.
* `<output_prefix>.frags.tsv.bgz` *(optional)*  
    - Fragment file corresponding to the downsampled BAM file.
* `<output_prefix>.summary.txt`  
    - Summary file describing the edit information and resulting statistics.
* `<output_prefix>.txt`  
    - List of selected read names.  
      Useful when storing a full BAM file is impractical. You can regenerate the BAM file later using this list:

```    
$ scBAMpler generateBAM \
    --input_bam test_data/HEPG2_subset.bam \
    --output_bam test_data/HEPG2_subset_c500_s12.bam \
    --selected_reads test_data/HEPG2_subset_c500_s12.txt

# ~XX min on subset (2.8Gb), ~XX min on full set (XX)
```

---------------

## Cell Homogeneity Extension

    $ command





---------------
## Citation

If you use scBAMpler in your research, please cite it!

- TBD

```
@article
```

