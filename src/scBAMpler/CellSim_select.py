"""
select-combos.py — Automatically Select Pseudo-bulk Combinations

Overview
--------
Takes the CSV output of `gen-pseudobulk-combos` and selects a representative
set of combinations by evenly sampling across two axes:
  - X axis: distance to a reference label centroid (chromatin similarity gradient)
  - Y axis: total peak read depth (sequencing depth bins)

This mirrors a manual inspection workflow but runs automatically, ensuring
even coverage of the parameter space rather than clustering selections in
low-distance / high-depth regions.

Optionally writes per-combo barcode lists for use with sinto filterbarcodes
(--write-barcodes).

Usage
-----
    # Basic: select combos and write summary stats
    python select-combos.py \\
        --input combos_all.csv \\
        --pickle medoids.pickle \\
        --output selected/ \\
        --ref-labels K562 HEPG2 \\
        --n-per-group 20

    # Also write per-combo barcode CSV files for sinto
    python select-combos.py \\
        --input combos_all.csv \\
        --pickle medoids.pickle \\
        --output selected/ \\
        --ref-labels K562 HEPG2 \\
        --n-per-group 20 \\
        --write-barcodes

Arguments
---------
    --input             CSV from gen-pseudobulk-combos (can be a concatenation
                        of multiple runs)
    --pickle            Pickle file from make-pseudobulks (needed to resolve
                        cluster → barcode mappings)
    --output            Output directory (created if it does not exist)
    --ref-labels        Reference labels to use as the X axis for selection.
                        One selection pass is run per label.
                        (default: all labels found in the data)
    --n-per-group       Number of combinations to select per reference label
                        per read-depth bin (default: 20)
    --read-depth-bins   Target read depths (Y axis bins) in peak read pairs.
                        (default: 5e7 1e8 1.5e8 2e8)
    --depth-tolerance   Half-window around each read depth bin; combos outside
                        [bin - tol, bin + tol] are excluded from that bin.
                        (default: 2e7)
    --max-sse           Exclude combinations with SSE above this threshold
                        (high SSE = poorly cohesive cluster). (default: 0.06)
    --max-intra-dist    For intra-label combinations (dominant_label_perc=100),
                        exclude if distance to ref label exceeds this value.
                        (default: 0.5)
    --label-col         Name of the grouping column (default: CellLine)
    --write-barcodes    If set, write one CSV per selected combo with columns
                        barcode and label, suitable for sinto filterbarcodes.

Output (in --output directory)
------
    selected_combos.csv     All selected combinations with full summary stats
    barcodes/               (if --write-barcodes) One CSV per combo:
                              combo_<ID>.barcodes.csv  with columns CB, <label_col>

Dependencies
------------
    pip install pandas numpy
"""

import os
import re
import ast
import pickle
import argparse

import numpy as np
import pandas as pd


# ── FILTERING ─────────────────────────────────────────────────────────────────

def apply_quality_filters(df, max_sse, max_intra_dist, ref_label):
    """
    Remove low-quality combinations before selection.

    Removes:
      - Combinations with SSE above max_sse (poorly cohesive)
      - Intra-label combinations (100% dominant) that are far from the
        reference label centroid (likely outlier clusters)

    Parameters
    ----------
    df : pd.DataFrame
    max_sse : float
    max_intra_dist : float
    ref_label : str

    Returns
    -------
    pd.DataFrame
    """
    dist_col = f"dist_to_{ref_label}"
    intra_mask = (
        (df["dominant_label_perc"] == 100.0) &
        (df[dist_col] > max_intra_dist)
    )
    sse_mask = df["sse_pearson_corr"] > max_sse
    removed = (intra_mask | sse_mask).sum()
    print(f"  [{ref_label}] Filtered out {removed} combos "
          f"(SSE>{max_sse}: {sse_mask.sum()}, intra+dist>{max_intra_dist}: {intra_mask.sum()})")
    return df[~(intra_mask | sse_mask)].copy()


# ── SELECTION ─────────────────────────────────────────────────────────────────

def bin_distance_axis(df, ref_label, n_bins=32):
    """
    Assign distance bins along the X axis (distance to ref label centroid).

    Uses denser bins at low distances where most variation occurs,
    sparser bins at high distances.

    Parameters
    ----------
    df : pd.DataFrame
    ref_label : str
    n_bins : int

    Returns
    -------
    pd.DataFrame with 'dist_bin' column added
    """
    dist_col = f"dist_to_{ref_label}"
    dist_max = df[dist_col].max()
    split = 0.35  # denser bins below this threshold

    bins_low = np.linspace(0, split, int(n_bins * 0.6))
    bins_high = np.linspace(split, dist_max, int(n_bins * 0.4) + 1)[1:]
    bins = np.concatenate([bins_low, bins_high])

    df = df.copy()
    df["dist_bin"] = np.digitize(df[dist_col], bins, right=True)
    return df


def select_points(df, ref_label, n_select, read_depth_bins, depth_tolerance):
    """
    Select n_select combinations per read depth bin by evenly sampling
    across the distance (X) axis.

    For each bin on the Y axis (read depth), picks the combination in each
    distance bin closest to that depth target. If a bin has no candidates
    within depth_tolerance, it is skipped. Random tie-breaking ensures
    diversity across runs.

    Parameters
    ----------
    df : pd.DataFrame
        Must have 'dist_bin' column from bin_distance_axis.
    ref_label : str
    n_select : int
        Combinations to select per read depth bin.
    read_depth_bins : list of float
    depth_tolerance : float

    Returns
    -------
    pd.DataFrame
        Selected rows with 'read_group' and 'ref_label' columns added.
    """
    dist_col = f"dist_to_{ref_label}"
    selected_rows = []

    for depth_target in read_depth_bins:
        lo = depth_target - depth_tolerance
        hi = depth_target + depth_tolerance

        candidates = df[df["total_peakreads"].between(lo, hi)].copy()
        if candidates.empty:
            print(f"  [{ref_label}] No candidates within ±{depth_tolerance:.0e} "
                  f"of depth={depth_target:.0e}, skipping")
            continue

        candidates["_dist_to_target"] = np.abs(
            candidates["total_peakreads"] - depth_target
        )

        # From each distance bin pick the candidate closest to depth target
        idx_best = (candidates
                    .groupby("dist_bin")["_dist_to_target"]
                    .idxmin()
                    .values)

        if len(idx_best) >= n_select:
            chosen = np.random.choice(idx_best, size=n_select, replace=False)
        else:
            chosen = idx_best

        sub = df.loc[chosen].copy()
        sub["read_group"] = f"{depth_target:.2e}"
        sub["ref_label"] = ref_label
        selected_rows.append(sub)
        print(f"  [{ref_label}] depth={depth_target:.0e}: selected {len(sub)}")

    if not selected_rows:
        return pd.DataFrame()
    return pd.concat(selected_rows, ignore_index=True)


# ── BARCODE OUTPUT ────────────────────────────────────────────────────────────

def write_barcodes(selected_df, embedding_df, label_col, out_dir):
    """
    Write one barcode CSV per selected combination for use with sinto.

    File format (no header):
        <barcode> <label>

    Parameters
    ----------
    selected_df : pd.DataFrame
    embedding_df : pd.DataFrame
    label_col : str
    out_dir : str
    """
    bc_dir = os.path.join(out_dir, "barcodes")
    os.makedirs(bc_dir, exist_ok=True)

    for idx, row in selected_df.iterrows():
        combo_id = row["ComboID"]
        medoids = row["group_medoids"]

        # Parse medoids from string representation if needed
        if isinstance(medoids, str):
            medoids = ast.literal_eval(medoids)

        cb_df = embedding_df[embedding_df["Cluster"].isin(medoids)][["CB", label_col]]

        # Strip any sample prefix (e.g. "SAMPLE#BARCODE" → "BARCODE")
        cb_df = cb_df.copy()
        cb_df["CB"] = cb_df["CB"].str.replace(r"^.*?#", "", regex=True)

        out_path = os.path.join(bc_dir, f"combo_{combo_id}.barcodes.csv")
        cb_df.to_csv(out_path, index=False, header=False, sep=" ")

    print(f"Wrote {len(selected_df)} barcode files to: {bc_dir}/")


# ── MAIN ──────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        prog="select-combos.py",
        description=(
            "Select a representative set of pseudo-bulk combinations by evenly "
            "sampling across distance-to-reference and read depth axes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, metavar="FILE",
        help="CSV from gen-pseudobulk-combos")
    parser.add_argument("--pickle", required=True, metavar="FILE",
        help="Pickle file from make-pseudobulks (for barcode resolution)")
    parser.add_argument("--output", required=True, metavar="DIR",
        help="Output directory")
    parser.add_argument("--ref-labels", nargs="+", default=None, metavar="LABEL",
        help="Reference labels for X-axis selection passes. "
             "Defaults to all labels found in dominant_label column.")
    parser.add_argument("--n-per-group", type=int, default=20, metavar="N",
        help="Combinations to select per reference label per read-depth bin")
    parser.add_argument("--read-depth-bins", type=float, nargs="+", metavar="N",
        default=[5e7, 1e8, 1.5e8, 2e8],
        help="Target read depths (Y axis bins)")
    parser.add_argument("--depth-tolerance", type=float, default=2e7, metavar="N",
        help="Half-window around each read depth bin")
    parser.add_argument("--max-sse", type=float, default=0.06, metavar="N",
        help="Exclude combinations with SSE above this value")
    parser.add_argument("--max-intra-dist", type=float, default=0.5, metavar="N",
        help="Exclude intra-label combos with distance to ref label above this value")
    parser.add_argument("--label-col", default="CellLine", metavar="COL",
        help="Name of the grouping column")
    parser.add_argument("--write-barcodes", action="store_true",
        help="Write per-combo barcode CSVs for sinto filterbarcodes")
    parser.add_argument("--seed", type=int, default=42,
        help="Random seed for tie-breaking")
    return parser.parse_args()


def main():
    args = parse_args()
    np.random.seed(args.seed)
    os.makedirs(args.output, exist_ok=True)

    # ── Load data ─────────────────────────────────────────────────────────────
    print(f"Loading combinations: {args.input}")
    combo_df = pd.read_csv(args.input)
    print(f"  {len(combo_df)} total combinations")

    print(f"Loading pickle: {args.pickle}")
    with open(args.pickle, "rb") as f:
        embedding_df, medoid_df, medoid_stats, cor_df = pickle.load(f)

    # ── Determine reference labels ────────────────────────────────────────────
    ref_labels = args.ref_labels or sorted(combo_df["dominant_label"].dropna().unique())
    print(f"Reference labels: {ref_labels}")

    # ── Selection pass per reference label ───────────────────────────────────
    all_selected = []
    combo_id = 0

    for ref_label in ref_labels:
        dist_col = f"dist_to_{ref_label}"
        if dist_col not in combo_df.columns:
            print(f"  WARNING: '{dist_col}' not found in input, skipping {ref_label}")
            continue

        print(f"\nProcessing reference label: {ref_label}")
        filtered = apply_quality_filters(
            combo_df, args.max_sse, args.max_intra_dist, ref_label
        )
        binned = bin_distance_axis(filtered, ref_label)
        selected = select_points(
            binned, ref_label, args.n_per_group,
            args.read_depth_bins, args.depth_tolerance
        )

        if selected.empty:
            print(f"  No combinations selected for {ref_label}")
            continue

        # Assign globally unique ComboIDs
        selected = selected.copy()
        selected["ComboID"] = range(combo_id, combo_id + len(selected))
        combo_id += len(selected)
        all_selected.append(selected)

    if not all_selected:
        print("No combinations were selected. Check your filters and input data.")
        return

    selected_df = pd.concat(all_selected, ignore_index=True)

    # ── Add cell line proportions ─────────────────────────────────────────────
    def cell_proportions(row):
        medoids = row["group_medoids"]
        if isinstance(medoids, str):
            medoids = ast.literal_eval(medoids)
        sub = embedding_df[embedding_df["Cluster"].isin(medoids)][args.label_col]
        return sub.value_counts().to_dict()

    print("\nComputing cell line proportions...")
    selected_df[f"{args.label_col}_proportions"] = selected_df.apply(
        cell_proportions, axis=1
    )

    # ── Save summary CSV ──────────────────────────────────────────────────────
    out_csv = os.path.join(args.output, "selected_combos.csv")
    selected_df.to_csv(out_csv, index=False)
    print(f"\nSaved {len(selected_df)} selected combinations to: {out_csv}")

    # ── Optionally write barcode files ────────────────────────────────────────
    if args.write_barcodes:
        print("Writing barcode files...")
        write_barcodes(selected_df, embedding_df, args.label_col, args.output)

    print("Done.")


if __name__ == "__main__":
    main()
