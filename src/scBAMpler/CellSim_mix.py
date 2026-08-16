"""
gen-pseudobulk-combos.py — Generate and Score Pseudo-bulk Combinations

Overview
--------
Takes the pickle output of `make-pseudobulks` and generates a large set of
candidate pseudo-bulk combinations by randomly sampling clusters, then scores
each combination on:
  - Total peak read depth
  - Mean pairwise Pearson correlation (chromatin similarity)
  - SSE of pairwise correlations (cluster cohesion)
  - Dominant label and its percentage
  - Distance to each reference label's centroid

Results are written to a CSV that can be passed to `select-combos` for
automated selection, or inspected manually in the companion notebook.

Usage
-----
    # Sample 2000 combinations from ALL clusters (unbiased)
    python gen-pseudobulk-combos.py \\
        --input medoids.pickle \\
        --output combos_all.csv \\
        --groups all \\
        --n-combos 2000 \\
        --cluster-size 500

    # Sample 1000 combinations from clusters dominated by K562 or HEPG2 only
    python gen-pseudobulk-combos.py \\
        --input medoids.pickle \\
        --output combos_k562_hepg2.csv \\
        --groups K562 HEPG2 \\
        --n-combos 1000 \\
        --cluster-size 500

    # Append multiple runs to build a mixed strategy
    cat combos_all.csv combos_k562_hepg2.csv > combos_combined.csv

Arguments
---------
    --input         Path to pickle file from make-pseudobulks
    --output        Path for output CSV file
    --groups        Labels to restrict sampling to (dominant label per cluster),
                    or 'all' for unbiased sampling across all clusters
    --n-combos      Number of random combinations to generate (default: 2000)
    --cluster-size  Cells per pseudo-bulk cluster, used to set combination size
                    range (default: 500)
    --ft-sizes      Target footprint sizes in cells to sample combinations for.
                    Combinations of r=ft_size/cluster_size clusters are drawn.
                    (default: 500 1000 2000 5000 10000 15000 20000)
    --label-col     Name of the grouping column (default: CellLine)
    --nproc         Number of parallel processes (default: 8)
    --seed          Random seed for reproducibility (default: 42)

Output
------
    CSV with one row per combination:
      group_medoids             list of cluster IDs in this combination
      total_peakreads           total raw read depth across clusters
      mean_pearson_corr         mean pairwise Pearson correlation
      sse_pearson_corr          SSE of pairwise correlations (cohesion)
      dominant_label            most common label in the combination
      dominant_label_perc       % of cells from the dominant label
      closest_label             label centroid closest to this combination
      label_dist_*              distance to each label's centroid (one col each)
      groups_sampled            value of --groups used to generate this row

"""

import os
import math
import pickle
import random
import argparse
import warnings
import functools
import multiprocessing as mp

import numpy as np
import pandas as pd

from datetime import datetime
from scipy.sparse import csc_matrix


# ── TIMER ─────────────────────────────────────────────────────────────────────

def timer(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = datetime.now()
        result = func(*args, **kwargs)
        print(f"--- {datetime.now() - start} h:m:s to run '{func.__name__}' ---")
        return result
    return wrapper


# ── COMBINATION GENERATION ────────────────────────────────────────────────────

def generate_combinations(medoid_arr, r, n, seed=42):
    """
    Draw n unique random combinations of size r from medoid_arr.

    Parameters
    ----------
    medoid_arr : list of str
        Pool of cluster IDs to sample from.
    r : int
        Number of clusters per combination.
    n : int
        Number of unique combinations to generate.
    seed : int
        Random seed.

    Returns
    -------
    list of tuple
    """
    if r > len(medoid_arr):
        return []
    random.seed(seed)
    combos = set()
    max_attempts = n * 20
    attempts = 0
    while len(combos) < n and attempts < max_attempts:
        combos.add(tuple(sorted(random.sample(medoid_arr, r))))
        attempts += 1
    return list(combos)


# ── COMBINATION SCORING ───────────────────────────────────────────────────────

def init_worker(medoid_stats, cor_df, embedding_df, label_col, all_labels):
    """Initialize shared state in each worker process."""
    global MEDOID_STATS, COR_DF, EMBEDDING_DF, LABEL_COL, ALL_LABELS
    MEDOID_STATS = medoid_stats
    COR_DF = cor_df
    EMBEDDING_DF = embedding_df
    LABEL_COL = label_col
    ALL_LABELS = all_labels


def process_combo(combo):
    """
    Score a single combination of clusters.

    Returns a list with:
      combo, total_peakreads, mean_pearson_corr, sse_pearson_corr,
      dominant_label, dominant_label_perc, closest_label, [dist_per_label...]
    """
    combo = np.array(combo)

    # Read depth
    depth = MEDOID_STATS.loc[combo, "PeakReads"].sum()

    # Pairwise correlations
    sub = COR_DF.loc[combo, combo].values
    corr_vect = sub[np.tril_indices(len(combo), k=-1)]
    mean_cor = np.mean(corr_vect)
    sse = np.sum(np.square(corr_vect - mean_cor)) / len(corr_vect)

    # Dominant label
    sub_labels = EMBEDDING_DF.loc[EMBEDDING_DF.Cluster.isin(combo), LABEL_COL]
    tmp = sub_labels.value_counts(normalize=True)
    dominant = tmp.idxmax()
    dom_perc = tmp.iloc[0] * 100

    # Distance to each label centroid (1 - mean correlation)
    dist_series = 1 - COR_DF.loc[ALL_LABELS, combo].mean(axis=1)
    closest_label = dist_series.idxmin()
    cl_dists = dist_series.values

    # return the cluster IDs as a plain list of str, not the ndarray. str(ndarray) is
    return [[str(c) for c in combo], depth, mean_cor, sse,
            dominant, dom_perc, closest_label, cl_dists]


@timer
def gather_combo_metrics(possible_combos, medoid_stats, cor_df, embedding_df,
                         label_col, all_labels, nproc=8):
    """
    Score all combinations in parallel.

    Returns
    -------
    pd.DataFrame
    """
    with mp.Pool(
        processes=nproc,
        initializer=init_worker,
        initargs=(medoid_stats, cor_df, embedding_df, label_col, all_labels)
    ) as pool:
        results = pool.map(process_combo, possible_combos)

    dist_cols = [f"dist_to_{l}" for l in all_labels]
    rows = []
    for r in results:
        row = r[:7] + list(r[7])
        rows.append(row)

    cols = ["group_medoids", "total_peakreads", "mean_pearson_corr",
            "sse_pearson_corr", "dominant_label", "dominant_label_perc",
            "closest_label"] + dist_cols
    return pd.DataFrame(rows, columns=cols)


# ── MAIN ──────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        prog="gen-pseudobulk-combos.py",
        description=(
            "Generate and score random pseudo-bulk combinations from cluster embeddings. "
            "Run multiple times with different --groups to build a mixed sampling strategy, "
            "then concatenate outputs."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, metavar="FILE",
        help="Pickle file from make-pseudobulks")
    parser.add_argument("--output", required=True, metavar="FILE",
        help="Output CSV file path")
    parser.add_argument("--groups", nargs="+", default=["all"], metavar="LABEL",
        help="Labels to restrict cluster pool to (dominant label per cluster), "
             "or 'all' for unbiased sampling. E.g.: --groups K562 HEPG2")
    parser.add_argument("--n-combos", type=int, default=2000, metavar="N",
        help="Number of random combinations to generate per ft-size")
    parser.add_argument("--cluster-size", type=int, default=500, metavar="N",
        help="Cells per cluster (from make-pseudobulks), used to compute r=ft_size/cluster_size")
    parser.add_argument("--ft-sizes", type=int, nargs="+", metavar="N",
        default=[500, 1000, 2000, 5000, 10000, 15000, 20000],
        help="Target footprint sizes in cells. One round of sampling per size.")
    parser.add_argument("--label-col", default="CellLine", metavar="COL",
        help="Name of the grouping column in embedding_df")
    parser.add_argument("--nproc", type=int, default=8, metavar="N",
        help="Number of parallel scoring processes")
    parser.add_argument("--seed", type=int, default=42,
        help="Random seed for reproducibility")
    return parser.parse_args()

# ── MAIN ──────────────────────────────────────────────────────────────────────


def main(args):
    #args = parse_args()

    #create the output directory if it doesn't exist, so steps can be run out of order
    out_dir = os.path.dirname(args.output)
    if out_dir: os.makedirs(out_dir, exist_ok=True)

    # ── Load pickle ───────────────────────────────────────────────────────────
    print(f"Loading: {args.input}")
    with open(args.input, "rb") as f:
        embedding_df, medoid_df, medoid_stats, cor_df = pickle.load(f)

    all_labels = sorted(embedding_df[args.label_col].unique())
    print(f"Labels found: {all_labels}")

    # ── Build cluster pool ────────────────────────────────────────────────────
    # medoid_stats holds one row per pseudo-bulk cluster AND one per label centroid.
    # Identify the centroids by name rather than by an "m" prefix: a label such as
    # "myeloid" would otherwise be mistaken for a cluster.
    label_set = set(all_labels)
    all_clusters = [i for i in medoid_stats.index if i not in label_set]

    # Assign dominant label to each cluster if not already present
    if "DominantLabel" not in medoid_stats.columns:
        dominant_labels = []
        for c in medoid_stats.index:
            if c not in label_set:
                sub = embedding_df.loc[embedding_df.Cluster == c, args.label_col]
                dominant_labels.append(sub.value_counts().idxmax() if len(sub) else None)
            else:
                dominant_labels.append(c)  # label centroids are their own dominant
        medoid_stats["DominantLabel"] = dominant_labels

    # Filter pool by --groups
    use_all = args.groups == ["all"]
    if use_all:
        cluster_pool = all_clusters
        group_tag = "all"
        print(f"Sampling from all {len(cluster_pool)} clusters")
    else:
        cluster_pool = medoid_stats[
            medoid_stats["DominantLabel"].isin(args.groups)
        ].index.tolist()
        cluster_pool = [c for c in cluster_pool if c not in label_set]
        group_tag = "+".join(sorted(args.groups))
        print(f"Sampling from {len(cluster_pool)} clusters dominated by: {args.groups}")

    if not cluster_pool:
        raise ValueError(
            f"No clusters found for groups: {args.groups}. "
            f"Available labels: {all_labels}"
        )

    # ── Generate combinations ─────────────────────────────────────────────────
    print(f"Generating combinations for ft-sizes: {args.ft_sizes}")

    #a combination of r clusters represents r * cluster_size cells, so the largest
    #population this dataset can express is the whole pool.
    pool_n = len(cluster_pool)
    max_ft = pool_n * args.cluster_size

    possible_combos = []
    for ft_size in args.ft_sizes:
        r = max(1, round(ft_size / args.cluster_size))

        if r > pool_n:
            msg = (
                f"ft-size {ft_size} needs {r} clusters but only {pool_n} are available, "
                f"so no combinations were generated for it. The largest population this "
                f"dataset supports is {max_ft} cells "
                f"({pool_n} clusters x {args.cluster_size} cells). Reduce --ft-sizes, "
                f"reduce --cluster-size in make-pseudobulks, or use a larger dataset."
            )
            warnings.warn(msg)
            print(f"  ft_size={ft_size:>6}  r={r}  generated=0   SKIPPED: {msg}")
            continue

        #with combinations sampled as sets, C(pool, r) is a hard ceiling on how many
        #distinct combinations can exist, regardless of --n-combos.
        n_possible = math.comb(pool_n, r)
        if n_possible < args.n_combos:
            warnings.warn(
                f"ft-size {ft_size} (r={r}) has only {n_possible} distinct combinations "
                f"available from {pool_n} clusters, fewer than the {args.n_combos} "
                f"requested. All {n_possible} will be returned."
            )

        combos = generate_combinations(cluster_pool, r=r, n=args.n_combos, seed=args.seed)
        possible_combos.extend(combos)
        print(f"  ft_size={ft_size:>6}  r={r}  generated={len(combos)}"
              f"{'' if n_possible >= args.n_combos else f'  (max possible: {n_possible})'}")

    if not possible_combos:
        raise ValueError(
            f"No combinations could be generated for any requested ft-size "
            f"{args.ft_sizes}. With {pool_n} clusters of {args.cluster_size} cells, "
            f"ft-sizes must be at most {max_ft}."
        )

    print(f"Total combinations to score: {len(possible_combos)}")

    # ── Score combinations ────────────────────────────────────────────────────
    combo_df = gather_combo_metrics(
        possible_combos, medoid_stats, cor_df, embedding_df,
        label_col=args.label_col, all_labels=all_labels, nproc=args.nproc
    )
    combo_df["groups_sampled"] = group_tag

    # ── Save ──────────────────────────────────────────────────────────────────
    combo_df.to_csv(args.output, index=False)
    print(f"Saved {len(combo_df)} combinations to: {args.output}")
    print("Done.")

