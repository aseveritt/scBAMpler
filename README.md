
<p align ="center">
<img src='https://github.com/aseveritt/scBAMpler/blob/main/docs/scBAMpler.png' style="max-width: 100%; height: auto;">
</p>

scBAMpler was developed to alter one aspect of a scATAC-seq dataset’s at a time: read count, cell count, and fraction of reads in peaks (FRiP) while preserving the original cell attributes. An extension was further developed which, using multiple cell types, will alter the cell-to-cell homogeneity of a population. 

![DOI](TBD)

## Installation

Please clone the repository:

    $ git clone git+https://github.com/aseveritt/scBAMpler.git

Then, create an environment with required dependencies. Installation and information about miniforge can be found [here](https://github.com/conda-forge/miniforge)

    $ miniforge create -n scBAMpler_env python=3.10 numpy scipy pandas jupyter samtools bamTools macs3
    $ pip install --editable .


---------------

## Download Test Data
```
$ cd test_data
$ scp ### HEPG2_subset.bam

```

---------------

## Data Quality Usage

### 1. Call Peak Locations

First, prepare your peak file. This can be done using any method you prefer—the only strict requirement is that the file be in **BED6 format**.
If you would like to use the peak standardization code from our manuscript, we provide it here:

```
$ bash helper_scripts/peak_calling/setup.sh
# ~20 min

$ Rscript helper_scripts/peak_calling/call_peaks.R \
    --bam_file test_data/HEPG2_subset.bam \
    --outdir test_data/ \
    --peak_length 500 \
    --cores 4

# ~XX min on subset (2.8Gb), ~XX min on full set (XX)
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
    --output_file test_data/HEPG2_subset.pickle

# ~10 min on subset (2.8Gb), ~XX min on full set (XX)
```
   
#### Input parameters  
* `--bam`
    - Path to the coordinate-sorted input BAM file.
* `--peak_file`
    - Path to the peak file in BED6 format.
* `--out_file`
    - Path where the final dictionary will be saved (as a `.pickle` file).

#### Output    
* `<output_file>`  
    - A Python pickle file containing a dictionary of all cell barcodes, their mapping to peak and non-peak reads, and the necessary numeric encoders. (e.g. HEPG2_subset.pickle)
* `<outfile>`.summary.txt
    - A plain-text file with summary statistics about the cell type.  (e.g. HEPG2_subset.summary.txt)
* `<outfile>`.reads_in_peaks.bed.gz  *(optionally, deleted with --delete_intersect flag)*
    - A gzipped BED-like file with two columns: cell barcode and associated peak-read QNAMEs.  
      *(Optionally deleted using the `--delete_intersect` flag.)*
    

### 3. Strategically Downsample BAM
Next, specify what feature to downsample and the extent. If you request a number larger than the cell population allows, it will give an error. 

```
$ scBAMpler sampler \
    --bam test_data/HEPG2_subset.bam \
    --peak_file test_data/? \
    --output_file test_data/HEPG2_subset.pickle \
```

#### Input parameters  
* --input_file
    - info
* --output_file
    - info
* --edit
    - info
* --value
    - info
* --seed
    - info
* --nproc
    - info
* --bam_file ??
    - info

#### Output 


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

