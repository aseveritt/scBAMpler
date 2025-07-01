# scBAMpler

<p align ="center">
<img src='https://github.com/aseveritt/scBAMpler/docs/scBAMpler.png' height="150">
</p>

[![DOI](TBD)
[![Docs](TBD)

scBAMpler was developed to alter one aspect of a scATAC-seq dataset’s at a time: read count, cell count, fraction of reads in peaks (FRiP), and cell-to-cell omogeneity while preserving the original cell attributes. 

## Installation

Please clone the repository:

    $ git clone git+https://github.com/aseveritt/scBAMpler.git

Then, create an environment with required dependencies. Installation and information about miniforge can be found [![here](https://github.com/conda-forge/miniforge)

    $ miniforge create -n scBAMpler_env python=3.10 numpy scipy pandas jupyter samtools bamTools macs3


---------------

## Data Quality Usage:

First, prepare your peak file. It is only strictly neccessary for the file to have XYZ. 
If you would like to use our mansucripts peak standardization code, we provide the code here, but it is not strictly neccessary. 
```
$ bash helper_scripts/peak_calling/setup.sh 
$ Rscript helper_scripts/peak_calling/call_peaks.R \
    --bam test_data/HEPG2_subset.bam \
    --outdir test_data/ \
    --peak_len 500 \
    --cores 4
```
##### Input parameters
* --bam
    - Path to BAM file
* --outdir
    - Path to output directory
* --peak_len
    - Standardize all peaks to this length
* --cores
    - optional, Number of cores available to mclapply
* --txdb
    - Txdb package used in get chromosomes lengths, default TxDb.Hsapiens.UCSC.hg38.knownGene
* --exclusion_list
    - List of regions to exclude from peak regions. 
* --summitfile
    - If MACS3 file already exists, just run standardization section. 

##### Output
* --bam
    - info


Next, build the dictionaries for each cell-type you would like to downsample:

```
$ scBAMpler create-dictionary \
    --bam test_data/HEPG2_subset.bam \
    --peak_file test_data/? \
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

## Cell Homogeneity Usage:

    $ command





---------------
## Citation

If you use scBAMpler in your research, please cite it!

- TBD

```
@article
```

