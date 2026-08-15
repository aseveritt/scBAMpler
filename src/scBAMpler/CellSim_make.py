"""
Steps
-----
1. Load HDF5 file (peak matrix + embedding)
2. Cluster cells within each cell line using size-constrained k-means
3. Compute centroid coordinates for each cluster and cell line
4. Normalize the peak matrix with TF-IDF
5. Average normalized accessibility within each cluster → pseudo-bulk profiles
6. Compute pairwise correlations across pseudo-bulk profiles
7. Save all results to a pickle file

The output is a pickle file containing:
  - embedding_df  : per-cell embedding coordinates with cluster assignments
  - medoid_df     : peak-by-pseudobulk accessibility matrix (TF-IDF normalized)
  - medoid_stats  : centroid coordinates and total read depth per pseudobulk
  - cor_df        : pairwise Pearson correlation matrix across pseudobulks

Arguments
---------
    --input         Path to input HDF5 file (see H5 structure in README)
    --output        Path for output pickle file
    --dimred        Embedding to use for clustering: 'umap' or 'tsne' (default: umap)
    --label-col     Name of the grouping column in the embedding (default: CellLine)
    --cluster-size  Target number of cells per pseudo-bulk cluster (default: 5000)
    --nproc         Number of parallel processes for clustering (default: 8)

Dependencies
------------
    pip install h5py pandas numpy scipy k-means-constrained
"""

import h5py
import pickle
import time
import functools
import argparse
import warnings

import pandas as pd
import numpy as np

from datetime import datetime
from scipy.sparse import csc_matrix, diags
from k_means_constrained import KMeansConstrained
from scipy.spatial.distance import cdist
from multiprocessing.dummy import Pool as ThreadPool


def timer(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = datetime.now()
        result = func(*args, **kwargs)
        print(f"--- {datetime.now() - start} h:m:s to run '{func.__name__}' ---")
        return result
    return wrapper


# ── DATA LOADING ──────────────────────────────────────────────────────────────

@timer
def loadH5(input_h5, dimred_method="tsne", label_col="CellLine"):
    """
    Load peak matrix and cell embedding from an HDF5 file.

    Parameters
    ----------
    input_h5 : str
        Path to HDF5 file (see README for expected structure).
    dimred_method : str
        Which embedding to load: 'umap' or 'tsne'.
    label_col : str
        Name to give the grouping column (column 3 of the embedding table).
        Defaults to 'CellLine'; set via --label-col to match your dataset.

    Returns
    -------
    input_df : pd.DataFrame
        Per-cell embedding with columns: embedding1, embedding2, CB, <label_col>.
    sparse_peak : scipy.sparse.csc_matrix
        Sparse peak-by-cell accessibility matrix.
    sparse_peak_colnames : np.ndarray of str
        Cell barcodes corresponding to columns of sparse_peak.
    """
    with h5py.File(input_h5, "r") as f:

        dimred_data = f["embedding/%s_df" % dimred_method][:]
        input_df = pd.DataFrame.from_records(dimred_data)
        input_df.columns = ["embedding1", "embedding2", "CB", label_col] #rename here for consistency
        input_df["CB"] = input_df["CB"].astype(str)
        input_df[label_col] = input_df[label_col].astype(str)

        # Load sparse peak matrix (CSC format)
        x = f["peak_matrix/x"][:].astype(np.int16)
        i = f["peak_matrix/i"][:]
        p = f["peak_matrix/p"][:]
        sparse_peak = csc_matrix((x, i, p))
        sparse_peak_colnames = f["peak_matrix/colnames"][:].astype(str)

    return input_df, sparse_peak, sparse_peak_colnames


# ── CLUSTERING ────────────────────────────────────────────────────────────────

def _cluster_one_cellline(df, cluster_size, seed, cellline):
    """
    Run size-constrained k-means on a single cell group's embedding.

    Clusters that exceed cluster_size after fitting have their outermost
    cells (farthest from centroid) reassigned to cluster -1 (unassigned).

    Parameters
    ----------
    df : pd.DataFrame
        Subset of the full embedding DataFrame for one cell line.
    cluster_size : int
        Target (minimum) cells per cluster.
    seed : int
        Random seed for reproducibility.
    cellline : str
        Cell line label (used only for logging).

    Returns
    -------
    df : pd.DataFrame or None
        DataFrame with a 'Cluster' column added, or None if too few cells.
    """
    
    X = df[["embedding1", "embedding2"]].to_numpy()

    if len(X) < cluster_size:
        warnings.warn(
            f"Skipping clustering for group '{cellline}': "
            f"only {len(X)} cells are available, which is fewer than the "
            f"requested cluster size of {cluster_size}. "
            f"All cells in this group will be marked as unassigned (Cluster = -1)."
        )
        df = df.copy()
        df["Cluster"] = -1
        return df

    # Constrained k-means clustering
    kmeans = KMeansConstrained(
        n_clusters=len(X) // cluster_size,
        size_min=cluster_size,
        size_max=None,
        init="k-means++",
        n_init=10,
        max_iter=300,
        random_state=seed,
        n_jobs=10,
    ).fit(X)

     # Assign cluster labels
    df["Cluster"] = kmeans.labels_
    cluster_centers = kmeans.cluster_centers_

    # Trim any clusters that still exceed the size limit
    cluster_counts = df["Cluster"].value_counts()
    large_clusters = cluster_counts[cluster_counts > cluster_size].index

    for c in large_clusters:
        mask = df["Cluster"] == c
        sub_df = df.loc[mask, ["embedding1", "embedding2"]]
        point_distances = cdist(
            sub_df.to_numpy(), [cluster_centers[c]], metric="euclidean"
        ).flatten()
        
        num_to_remove = len(point_distances) - cluster_size
        removal_indices = np.argpartition(point_distances, -num_to_remove)[-num_to_remove:]
        idxs_large = [sub_df.index[i] for i in removal_indices]
        df.loc[idxs_large, "Cluster"] = -1

    return df


@timer
def constrained_cluster(df, cluster_size=500, seed=42, nproc=8, label_col="CellLine"):
    """
    Cluster cells within each group using size-constrained k-means,
    run in parallel across groups.

    Cells that cannot be assigned to a full cluster are labeled 'm-1'.
    All other clusters receive a globally unique label (m1, m2, ...).

    Parameters
    ----------
    df : pd.DataFrame
        Full embedding DataFrame with columns: embedding1, embedding2, CB, <label_col>.
    cluster_size : int
        Target cells per cluster (default: 500).
    seed : int
        Random seed (default: 42).
    nproc : int
        Number of parallel threads (default: 8).
    label_col : str
        Name of the grouping column (default: 'CellLine').

    Returns
    -------
    df : pd.DataFrame
        Input DataFrame with 'Cluster' column added.
    """
    
    grouped = df.groupby(label_col)

    with ThreadPool(processes=nproc) as pool:
        results = pool.starmap(
            _cluster_one_cellline,
            [(df_sub, cluster_size, seed, cl) for cl, df_sub in grouped],
        )

    df = pd.concat([r for r in results if r is not None], ignore_index=True)
    df["Cluster"] = df["Cluster"].astype(str)

    # Assign globally unique cluster labels
    i = 1
    for (label, sub_df) in df.groupby([label_col, "Cluster"]):
        if label[1] == "-1":
            df.loc[sub_df.index, "Cluster"] = "m-1"
        else:
            df.loc[sub_df.index, "Cluster"] = f"m{i}"
            i += 1

    return df


# ── EMBEDDING COORDINATES ─────────────────────────────────────────────────────

@timer
def get_embedding_coord(embedding_df, label_col="CellLine"):
    """
    Compute centroid coordinates in embedding space for each group
    and each cluster.

    Parameters
    ----------
    embedding_df : pd.DataFrame
        Embedding DataFrame with columns: embedding1, embedding2, <label_col>, Cluster.
    label_col : str
        Name of the grouping column (default: 'CellLine').

    Returns
    -------
    info_df : pd.DataFrame
        One row per group / cluster with columns: Cluster, Kcenter_x, Kcenter_y.
    """
    
    info = {}

    for cl, df_sub in embedding_df.groupby(label_col):
        cx, cy = df_sub[["embedding1", "embedding2"]].mean().values
        info[cl] = [cl, cx, cy]

    for cl, df_sub in embedding_df.groupby("Cluster"):
        if cl == "m-1":
            continue
        cx, cy = df_sub[["embedding1", "embedding2"]].mean().values
        info[cl] = [cl, cx, cy]

    info_df = pd.DataFrame.from_dict(info, orient="index")
    info_df.columns = ["Cluster", "Kcenter_x", "Kcenter_y"]
    return info_df


# ── TF-IDF NORMALIZATION ──────────────────────────────────────────────────────

@timer
def perform_tfidf(csr_mat):
    """
    Apply TF-IDF normalization to a sparse peak-by-cell matrix.

    Follows the approach from Hill (2019) for scATAC-seq dimensionality
    reduction (https://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/).

    Term Frequency (TF): column-wise normalization (per cell), then log-scaled
    to approximate log-CPM. Inverse Document Frequency (IDF): row-wise, penalizes
    peaks open in many cells.

    Parameters
    ----------
    csr_mat : scipy.sparse matrix
        Raw peak-by-cell count matrix.

    Returns
    -------
    scipy.sparse.csr_matrix
        TF-IDF normalized matrix, same shape as input.
    """
    csr_mat = csr_mat.tocsr()
    #let csr_mat.shape = (10,300)
    
    # Term Frequency (TF)
    col_sums = np.array(csr_mat.sum(axis=0)).flatten()  # (1, 3000) -> (3000,)
    col_sums[col_sums == 0] = 1  # avoid division by zero
    tf = csr_mat.multiply(1 / col_sums) 
    #(10, 300)
    
    #Log transform
    tf.data = np.log1p(tf.data * 100000) #something at something similar to a log cpm now 
    #(10, 300)
    
    ## Inverse Document Frequency (IDF)
    row_sums = np.array(csr_mat.sum(axis=1)).flatten()  
    row_sums[row_sums == 0] = 1  # avoid division by zero
    idf = np.log1p(csr_mat.shape[0] / row_sums)  # Compute IDF
    idf = diags(idf)  # Convert to a sparse diagonal matrix
    #(10, 1)
        
    # Final TF-IDF transformation
    tfidf = idf @ tf
    #(10, 300) = (10, 300) * (10, 1)
    
    return tfidf.tocsr()


# ── PSEUDO-BULK PROFILES ──────────────────────────────────────────────────────

@timer
def get_peakmat_medoid(sparse_peak_norm, sparse_peak, sparse_peak_colnames, embedding_df, label_col="CellLine"):
    """
    Compute pseudo-bulk accessibility profiles by averaging normalized peak
    accessibility across all cells in each cluster (and each group).

    Parameters
    ----------
    sparse_peak_norm : scipy.sparse matrix
        TF-IDF normalized peak-by-cell matrix.
    sparse_peak : scipy.sparse matrix
        Raw peak-by-cell matrix (used for read depth calculation).
    sparse_peak_colnames : np.ndarray of str
        Cell barcodes corresponding to columns of sparse_peak.
    embedding_df : pd.DataFrame
        Embedding DataFrame with <label_col> and Cluster columns.
    label_col : str
        Name of the grouping column (default: 'CellLine').

    Returns
    -------
    medoid_df : pd.DataFrame
        peaks × pseudobulks matrix of averaged TF-IDF accessibility.
    read_depth : dict
        Total raw read count per pseudobulk label.
    """
    
    medoid_dict = {}
    read_depth = {}
    colnames_to_index = {name: idx for idx, name in enumerate(sparse_peak_colnames)}

    embedding_df_grouped = embedding_df.groupby(label_col)
    for c, df_sub in embedding_df_grouped:
        keep_cbs = df_sub["CB"].values #get all cell barcodes belonging to that cell line. 
        idxs = [colnames_to_index[cb] for cb in keep_cbs] # if cb in colnames_to_index]

        medoid_dict[c] = np.asarray(sparse_peak_norm[:, idxs].mean(axis=1)).flatten()
        read_depth[c] = sparse_peak[:, idxs].sum() 
    
    embedding_df_grouped = embedding_df.groupby("Cluster")
    for c, df_sub in embedding_df_grouped:
        if c == "m-1": continue
        keep_cbs = df_sub["CB"].values #get all cell barcodes belonging to that cell line. 
        idxs = [colnames_to_index[cb] for cb in keep_cbs] # if cb in colnames_to_index]
        
        medoid_dict[c] = np.asarray(sparse_peak_norm[:, idxs].mean(axis=1)).flatten()
        read_depth[c] = sparse_peak[:, idxs].sum() 
    
    medoid_df = pd.DataFrame.from_dict(medoid_dict, orient='index')
    
    return medoid_df.T, read_depth


# ── MAIN ──────────────────────────────────────────────────────────────────────


def main(args):
    #args = parse_args()

    print(f"Loading HDF5: {args.input}")
    embedding_df, sparse_peak, sparse_peak_colnames = loadH5(
        args.input, args.dimred, label_col=args.label_col
    )

    print(f"Clustering cells (cluster size={args.cluster_size}, nproc={args.nproc})")
    embedding_df = constrained_cluster(
        embedding_df, cluster_size=args.cluster_size, nproc=args.nproc,
        label_col=args.label_col
    )

    print("Computing embedding centroids")
    medoid_stats = get_embedding_coord(embedding_df, label_col=args.label_col)

    print("Applying TF-IDF normalization")
    sparse_peak_norm = perform_tfidf(sparse_peak)

    print("Generating pseudo-bulk profiles")
    medoid_df, read_depth = get_peakmat_medoid(
        sparse_peak_norm, sparse_peak, sparse_peak_colnames, embedding_df,
        label_col=args.label_col
    )
    medoid_stats["PeakReads"] = read_depth

    print("Computing pairwise correlations")
    corr_matrix = np.corrcoef(medoid_df.T)
    cor_df = pd.DataFrame(corr_matrix, index=medoid_stats.index, columns=medoid_stats.index)

    print(f"Saving results to: {args.output}")
    with open(args.output, "wb") as f:
        pickle.dump([embedding_df, medoid_df, medoid_stats, cor_df], f)

    print("Done.")

