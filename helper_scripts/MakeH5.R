# makeH5.R
# Purpose: Extract peak matrix and embeddings from an ArchR project and save them to an HDF5 (.h5) file for downstream use (e.g. in Python).
#
# Output: An HDF5 file with the following structure (see peakmat_input_H5_README.md)
#   /peak_matrix/x         - non-zero values of the sparse peak matrix
#   /peak_matrix/i         - row indices (CSC format)
#   /peak_matrix/p         - column pointers (CSC format)
#   /peak_matrix/colnames  - cell barcodes (column names)
#   /embedding/umap_df     - UMAP coordinates + cell barcode + cell line
#   /embedding/tsne_df     - tSNE coordinates + cell barcode + cell line

# ── SETUP ─────────────────────────────────────────────────────────────────────
 
library(ArchR)
library(ggrepel)
library(dplyr)
library(pals)
library(rhdf5)
library(Matrix)
 
set.seed(1)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

peak_bed_file <- "test_data/union_standardized_500bp.bed"

h5_output_file <- "test_data/peakmat_input.h5"


# ── LOAD PROJECT AND PEAKS ────────────────────────────────────────────────────

arrows <- createArrowFiles(
  inputFiles = c("test_data/HEPG2_subset.bam", "test_data/K562_subset.bam", "test_data/MCF7_subset.bam"),
  sampleNames = c("HEPG2", "K562", "MCF7"),
  addTileMat = TRUE,
  addGeneScoreMat = FALSE,
  bcTag = "CB", 
  force=TRUE,
  minTSS = 0,
  minFrags = 0
)

proj <- ArchRProject(
  ArrowFiles = arrows,
  outputDirectory = "test_data/CLsub",
  copyArrows = TRUE
)

peak.gr <- import(peak_bed_file)
validChr <- as.character(seqnames(getChromSizes(proj)))
peak.gr <- peak.gr[seqnames(peak.gr) %in% validChr]
proj <- addPeakSet(ArchRProj = proj, peakSet = peak.gr, force = TRUE)
proj <- addPeakMatrix(proj)

proj$barcodes <- gsub(".*#","",proj$cellNames)
proj$CellLine <- sub("_.*", "", proj$Sample)

proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 4,
    clusterParams = list( #See Seurat::FindClusters
        resolution = 0.1,
        sampleCells = 10000,
        n.start = 10
    ),
    varFeatures = 15000,
    dimsToUse = 1:30
)

proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.1
)

proj <- addUMAP(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
)

proj <- addTSNE(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "TSNE",
    perplexity = 30
)


# ── EXTRACT PEAK MATRIX ───────────────────────────────────────────────────────
peakmat <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
sparse_peak_mat <- assay(peakmat)   # peaks x cells, stored in CSC (dgCMatrix) format
 
# ── WRITE HDF5 ────────────────────────────────────────────────────────────────
# Remove existing file if present, then create fresh
if (file.exists(h5_output_file)) file.remove(h5_output_file)
h5createFile(h5_output_file)
 
# --- peak_matrix group: sparse matrix in CSC format ---
# Reconstruct in Python/scipy with:
#   scipy.sparse.csc_matrix((x, i, p), shape=(n_peaks, n_cells))
h5createGroup(h5_output_file, "peak_matrix")
h5write(sparse_peak_mat@x, h5_output_file, "peak_matrix/x")   # non-zero values
h5write(sparse_peak_mat@i, h5_output_file, "peak_matrix/i")   # row indices
h5write(sparse_peak_mat@p, h5_output_file, "peak_matrix/p")   # column pointers
h5write(colnames(sparse_peak_mat), h5_output_file, "peak_matrix/colnames")  # cell barcodes
 
# --- embedding group: dimensionality reduction results ---
h5createGroup(h5_output_file, "embedding")
 
# UMAP: columns are IterativeLSI#UMAP_Dimension_1, IterativeLSI#UMAP_Dimension_2, CB, CellLine
tmp <- proj@embeddings$UMAP$df
tmp$CB <- rownames(tmp)
tmp$CellLine <- as.vector(proj@cellColData$CellLine)
h5write(tmp, h5_output_file, "embedding/umap_df")
 
# tSNE: same structure as UMAP
tmp <- proj@embeddings$TSNE$df
tmp$CB <- rownames(tmp)
tmp$CellLine <- as.vector(proj@cellColData$CellLine)
h5write(tmp, h5_output_file, "embedding/tsne_df")
 
# Close all HDF5 connections
h5closeAll()
 
message("Done! HDF5 written to: ", h5_output_file)
