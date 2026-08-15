
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
    $ pip install h5py k-means-constrained #specific to cell similarity extension
    $ cd scBAMpler/
    $ pip install .


---------------

## Download Test Data
Data is available on [Zenodo](TBD)
```
$ cd scBAMpler
$ wget https://zenodo.org/records/XXXXX/files/test_data.tar.gz?download=1
$ tar -xzvf test_data.tar.gz
#2.94Gb in size

# optionally, if you would like to see the example output directory without running the tutorial:
$ cd example_output/
$ wget https://zenodo.org/records/XXXXX/files/example_output.tar.gz?download=1
$ tar -xzvf example_output.tar.gz
#1.68Gb in size 
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
We assume the BAM file contains a cell barcode tag in the form `CB:Z:*` and botht the bam and peak file are coordinate sorted. 

```
$ scBAMpler create-dictionary \
    --bam test_data/HEPG2_subset.bam \
    --peak_file test_data/HEPG2_subset_standardized_500bp.bed \
    --output_file example_output/HEPG2_subset.pickle \
    --verbose
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


### 4. Recreate BAM from stored text file. 
It is unlikely users will want to store the subset bams in long-terms storage. One option is to store the .txt files only long term and recreate the downsampled bam files directly if needed again in the future.   
```    
$ scBAMpler generateBAM \
    --input_bam test_data/HEPG2_subset.bam \
    --output_bam example_output/HEPG2_subset_c500_s12_REMADE.bam \
    --selected_reads example_output/HEPG2_subset_c500_s12.txt \
    --nproc 5

$ diff <(samtools view example_output/HEPG2_subset_c500_s12.bam) <(samtools view example_output/HEPG2_subset_c500_s12_REMADE.bam) #returns nothing

# ~2 min on subset (2.8Gb), ~3 min on full set (25Gb)
```

---------------

## Cell Homogeneity Extension

### 1. Generate H5 input file
First, generate an HDF5 file to store the single-cell ATAC-seq data: a sparse peak-by-cell accessibility matrix plus UMAP and tSNE embeddings. The helper script creates this from an ArchR project, but any file following the structure below will work.

peakmat_input.h5
├── peak_matrix/
│   ├── x          float64[]   Non-zero values of the sparse matrix
│   ├── i          int32[]     Row indices (0-based) of non-zero values
│   ├── p          int32[]     Column pointers (length = n_cells + 1)
│   └── colnames   string[]    Cell barcodes, one per column (length = n_cells)
└── embedding/
    ├── umap_df    table       UMAP coordinates per cell (see below) OR
    └── tsne_df    table       tSNE coordinates per cell (see below)
    
Each embedding table has the form: (x coordinate, y coordinate, cell barcode, label-col)
To generate using the helper script:

```
$ Rscript helper_scripts/makeH5.R

```

### 2. Make small, pseudobulks of identical size. 
Next, we need to build pseudobulk profiles and collect summary information that will support the bottom-up approach for constructing mixed synthetic populations. 

```
$ scBAMpler make-pseudobulks \
    --input test_data/peakmat_input.h5 \
    --output example_output/medoids_s5000.pickle \
    --dimred umap \
    --label-col CellLine \
    --cluster-size 500 \
    --nproc 5

# ~XX min on subset (2.8Gb), ~XX min on full set (25Gb)
```
                               
#### Input Parameters
* `--input`  
    - Path to input HDF5 file (peakmat_input.h5)
* `--output`  
    - Path for output pickle file (e.g. medoids_s5000.pickle)
* `--dimred`  
    - Embedding to use for clustering 
    - Choices: `"umap"` or `"tsne"` (whatever is present in h5 file) 
* `--label-col`  
    - Name of the grouping column in the embedding (column 3 of the H5 embedding table)
* `--cluster-size`  
    - Target number of cells per pseudo-bulk cluster (default 500 cells)
* `--nproc`  
    - Number of parallel processes for clustering. 

#### Output
* `*.pickle (specified by --output)`  
    - Python object containing four items:
        1. `embedding_df` — per-cell embedding coordinates with cluster assignments
        2. `medoid_df` — peaks × pseudobulks matrix of TF-IDF normalized accessibility
        3. `medoid_stats` — centroid coordinates and total read depth per pseudobulk
        4. `cor_df` — pairwise Pearson correlation matrix across pseudobulks


### 3. Generate Hypothetical Cell Populations
DESCRIPTION TEXT .

```
#Sample 2000 combinations from all clusters
$ scBAMpler mix-pseudobulks \
    --input example_output/medoids.pickle \
    --output example_output/combos_all.csv \ 
    --groups all \
    --n-combos 2000 \
    --cluster-size 500

#Sample 1000 combinations from K562 and HEPG2 dominated clusters
$ scBAMpler mix-pseudobulks \
    --input example_output/medoids_s5000.pickle \
    --output example_output/combos_k562_hepg2.csv \
    --groups K562 HEPG2 \
    --n-combos 1000 \
    --cluster-size 500

cat combos_all.csv combos_k562_hepg2.csv > combos_combined.csv

# ~XX min on subset (2.8Gb), ~XX min on full set (25Gb)
```

#### Input Parameters
* `--input`  
    -  Path to pickle file from make-pseudobulks
* `--output`  
    - Path for output CSV file
* `--groups`  
    - Labels to restrict sampling to (dominant label per cluster), or 'all' for unbiased sampling across all clusters
    - Choices: `"all"` or comma separated list of label sets to specifically mix e.g. `"K562, HEPG2"`
* `--n-combos`  
    - Number of random combinations to generate (default: 2000)
* `--cluster-size`  
    - Cells per pseudo-bulk cluster
    - Should be a muliple of the clusters layed out above-- amanda I think we can infer this directly.
* `--ft-sizes`
    - Target footprint sizes in cells to sample combinations for. Combinations of r=ft_size/cluster_size clusters are drawn.
    - (default: 500 1000 2000 5000 10000 15000 20000)
* `--label-col`  
    - Name of the grouping column (default: CellLine)
* `--nproc`  
    - Number of parallel processes (default: 8)
* `--seed`  
    - Random seed for reproducibility (default: 42)

#### Output
* `*.csv (specified by --output)`
    - CSV with one row per hypothetical cell population:
      group_medoids             list of included cluster IDs
      total_peakreads           total raw read depth
      mean_pearson_corr         mean pairwise Pearson correlation
      sse_pearson_corr          SSE of pairwise correlations (cohesion)
      dominant_label            most common label
      dominant_label_perc       % of cells from the dominant label
      closest_label             label centroid closest to this combination
      label_dist_*              distance to each label's centroid (one col each)
      groups_sampled            value of --groups used to generate this row


### 4. Select Cell Populations of interest and optionally write cell-barcodes for downstream analyses. 
DESCRIPTION TEXT .
        
```
#Alternatively, if you want to visualize the distribution in a notebook and output from there you could run something similar to
helper_scripts/inspect-select-combos.ipynb

$ scBAMpler select-populations \
    --input example_output/combos_all.csv \
    --pickle example_output/medoids.pickle \ 
    --output selected/ \
    --ref-labels K562 HEPG2 \\ ???? AMANDA
    --n-per-group 20


# ~XX min on subset (2.8Gb), ~XX min on full set (25Gb)
```

#### Input Parameters
* `--input`  
    -  CSV from gen-pseudobulk-combos
* `--pickle`  
    - Pickle file from make-pseudobulks (needed to resolve cluster → barcode mappings)
* `--output`
    - Output directory (created if it does not exist)

* `--ref-labels`
    - Reference labels to use as the X axis for selection.
                        One selection pass is run per label.
                        (default: all labels found in the data)

* `--n-per-group`
    - Number of combinations to select per reference label per read-depth bin (default: 20)

* `--write-barcodes`
    - If set, write one CSV per selected combo with columns barcode and label, suitable for sinto filterbarcodes.

#### Output
* `selected_combos.csv`
    - All selected combinations with full summary stats
* `barcodes/`
    - (if --write-barcodes) One CSV per combo: combo_<ID>.barcodes.csv  with columns CB, <label_col>

      
---------------
## Citation

If you use scBAMpler in your research, please cite it!

```
@article{TBD}
```

