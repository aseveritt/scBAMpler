
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
# ~10 min on subset (2.8Gb), ~45 min on full set (14.Gb)
```

#### Input Parameters
* `--bam_file`  
    - Path to the input BAM file.
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
    - 
#### Input parameters
* --bam_file
    - Path to BAM file
* --outdir
    - Path to output directory
* --peak_length
    - Standardize all peaks to this length
* --cores
    - optional, Number of cores available to mclapply
* --txdb
    - Txdb package used in get chromosomes lengths, default TxDb.Hsapiens.UCSC.hg38.knownGene
* --exclusion_file
    - List of regions to exclude from peak regions. 
* --summit_file
    - If MACS3 file already exists, just run standardization section. 

#### Output
* /outdir/*_standardized_<peak_length>bp.bed
    - BED6 formatted output 


### 2. Build Cell Type Input Dictionaries
Next, build the dictionaries for each cell-type you would like to downsample:

```
$ scBAMpler create-dictionary \
    --bam test_data/HEPG2_subset.bam \
    --peak_file test_data/HEPG2_subset_standardized_500bp.bed \
    --output_file test_data/HEPG2_subset.pickle \
```
   
#### Input parameters  
* --bam
    - Path to BAM file
* --peak_file
    - Path to peak file in BED4 file? 
* --out_file
    - Name and location of final cell type dictionary stored as a pickle file
* --intersect_file
    - optional, provide an external file that maps reads to peaks. Must be in the form XYZ. 
* --delete_intersect
    - optional, by default the intersect file discussed above is stored. This flag will automatically delete it. 

#### Output    

In the output folder given in `--outdir`, the following files and folder structure will be created:
    - \<outdir\>/\<TF\>_overview.pickle
      This is a pickle file, or a compressed python file, that contains the dictionary of all cell barcodes, and the required encoders to numerically
    - \<outdir\>idk.stats
    

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

