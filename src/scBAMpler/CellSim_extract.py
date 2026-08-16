"""
CellSim_extract.py — Write Extraction Scripts for Selected Populations

Overview
--------
Takes a directory of barcode files describing selected populations and emits one bash
script per population. Each script filters every contributing label's BAM down to that
population's barcodes, merges the results into a single BAM, indexes it, and removes the
intermediates.

Nothing is executed here. Extraction is expensive -- one pass over a BAM per label per
population -- so the scripts are written for you to run, submit to a scheduler, or
inspect first.

A mixed population draws cells from more than one experiment, so its reads live in more
than one BAM. That is why one BAM must be given per label, and why each script runs sinto
once per label before merging.

Input
-----
There is no manifest file. Populations are discovered from the barcode filenames, which
carry everything needed:

    <barcode-dir>/combo_<ID>.<label>.barcodes.csv

  <ID>     integer identifying the population. One output BAM is written per ID.
  <label>  the group the cells belong to, matching a --bam LABEL.

Each file has one line per cell, two space separated columns and no header, which is what
`sinto filterbarcodes` expects for --cells:

    ATGATAGGACCTAGGC HEPG2
    AACTAGCACCGATCGC HEPG2

Barcodes must match the CB tag in the BAM, so any "<Sample>#" prefix has to be stripped
before writing them. notebooks/inspect-select-combos.ipynb writes this layout directly.

Arguments
---------
    --barcode-dir       Directory of barcode files, as described above
    --output            Output directory for scripts and BAMs
    --bam               LABEL=PATH for each label, repeated. e.g.
                        --bam HEPG2=test_data/HEPG2_subset.bam
                        --bam K562=test_data/K562_subset.bam
    --combo-ids         Only write scripts for these IDs (default: all found)
    --nproc             Processors passed to sinto and samtools (default: 8)
    --output_fragment   Also produce a .frags.tsv.bgz per population

Output (in --output)
------
    scripts/combo_<ID>.sh       One script per population, with absolute paths
    run_all.sh                  Runs every script above, in order
    combo_<ID>.bam (+ .bai)     Written when a script is run
    combo_<ID>.frags.tsv.bgz    (if --output_fragment)
"""

import os
import re
import sys
import glob
import stat
import argparse
import collections


BARCODE_RE = re.compile(r"^combo_(\d+)\.(.+)\.barcodes\.csv$")


def _fail(msg):
    print(f"ERROR: {msg}")
    sys.exit(1)


# ── INPUT DISCOVERY ───────────────────────────────────────────────────────────

def parse_bam_args(bam_args):
    """Turn ["HEPG2=a.bam", "K562=b.bam"] into {"HEPG2": "a.bam", ...}."""
    mapping = {}
    for item in bam_args:
        if "=" not in item:
            _fail(f"--bam expects LABEL=PATH, got '{item}'. "
                  f"e.g. --bam HEPG2=test_data/HEPG2_subset.bam")
        label, path = (p.strip() for p in item.split("=", 1))
        if not label or not path:
            _fail(f"--bam expects LABEL=PATH, got '{item}'")
        if not os.path.isfile(path):
            _fail(f"BAM for label '{label}' does not exist: {path}")
        mapping[label] = os.path.abspath(path)
    return mapping


def discover_populations(barcode_dir):
    """
    Read the barcode directory and return {combo_id: {label: barcode_file}}.

    The filenames are the source of truth for which labels a population draws from,
    so a label can never be silently omitted from an extraction.
    """
    found = collections.defaultdict(dict)
    ignored = []

    for path in sorted(glob.glob(os.path.join(barcode_dir, "*"))):
        name = os.path.basename(path)
        m = BARCODE_RE.match(name)
        if not m:
            if os.path.isfile(path):
                ignored.append(name)
            continue
        combo_id, label = int(m.group(1)), m.group(2)
        found[combo_id][label] = os.path.abspath(path)

    return dict(found), ignored


# ── SCRIPT GENERATION ─────────────────────────────────────────────────────────

def write_combo_script(script_path, combo_id, per_label, out_dir, nproc, output_fragment):
    """
    Write the bash script for one population.

    per_label : list of (label, bam_path, barcode_file), all absolute
    """
    tmp_dir = os.path.join(out_dir, "tmp", f"combo_{combo_id}")
    final_bam = os.path.join(out_dir, f"combo_{combo_id}.bam")

    L = ["#!/bin/bash",
         f"# Extract population combo_{combo_id} from {len(per_label)} label(s): "
         f"{', '.join(l for l, _, _ in per_label)}",
         "#",
         "# Requires sinto and samtools on PATH (conda activate scBAMpler_env).",
         "",
         "set -eo pipefail",
         "",
         f'TMP="{tmp_dir}"',
         f'FINAL="{final_bam}"',
         'mkdir -p "$TMP"',
         ""]

    # sinto names its output after the label in column 2 of the barcode file
    L.append("# ── filter each label's BAM to this population's barcodes ──")
    parts = []
    for label, bam, bc in per_label:
        L.append(f'echo "--> {label}"')
        L.append(f'sinto filterbarcodes -b "{bam}" -c "{bc}" -p {nproc} --outdir "$TMP"')
        parts.append(os.path.join(tmp_dir, f"{label}.bam"))
    L.append("")

    L.append("# ── merge into a single BAM ──")
    if len(parts) == 1:
        L.append(f'mv "{parts[0]}" "$FINAL"')
    else:
        L.append(f'samtools merge -f -o "$FINAL" --threads {nproc} '
                 + " ".join(f'"{p}"' for p in parts))
    L.append('samtools index "$FINAL"')
    L.append("")

    if output_fragment:
        frag = os.path.join(out_dir, f"combo_{combo_id}.frags.tsv.bgz")
        L.append("# ── fragment file ──")
        L.append(f'sinto fragments --collapse_within -p {nproc} -b "$FINAL" '
                 f'-f "$TMP/frags.tsv" > /dev/null')
        L.append(f'bedtools sort -i "$TMP/frags.tsv" '
                 f"""| awk '{{print $1, $2, $3, "combo_{combo_id}:"$4, $5}}' """
                 f"""| tr ' ' '\\t' | bgzip -c > "{frag}" """)
        L.append("")

    L += ["# ── clean up intermediates ──",
          'rm -rf "$TMP"',
          "",
          f'echo "--> done: {final_bam}"',
          ""]

    with open(script_path, "w") as f:
        f.write("\n".join(L))
    os.chmod(script_path, os.stat(script_path).st_mode | stat.S_IXUSR | stat.S_IXGRP)


def write_runner(runner_path, script_paths):
    """Write a script that runs every per-combo script in turn."""
    L = ["#!/bin/bash",
         f"# Runs all {len(script_paths)} population extraction scripts.",
         "# Each is independent, so these can also be submitted separately.",
         "",
         "set -eo pipefail",
         ""]
    for i, p in enumerate(script_paths, 1):
        L.append(f'echo "[{i}/{len(script_paths)}] {os.path.basename(p)}"')
        L.append(f'bash "{p}"')
    L += ["", f'echo "All {len(script_paths)} populations extracted."', ""]

    with open(runner_path, "w") as f:
        f.write("\n".join(L))
    os.chmod(runner_path, os.stat(runner_path).st_mode | stat.S_IXUSR | stat.S_IXGRP)


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main(args):

    bam_map = parse_bam_args(args.bam)
    print(f"BAMs provided for {len(bam_map)} label(s): {sorted(bam_map)}")

    barcode_dir = os.path.abspath(args.barcode_dir)
    if not os.path.isdir(barcode_dir):
        _fail(f"--barcode-dir is not a directory: {barcode_dir}")

    populations, ignored = discover_populations(barcode_dir)
    if not populations:
        _fail(f"No barcode files matching 'combo_<ID>.<label>.barcodes.csv' in {barcode_dir}. "
              f"See the module docstring for the expected layout."
              + (f" Ignored {len(ignored)} other file(s), e.g. {ignored[0]}." if ignored else ""))

    all_labels = sorted({lab for labs in populations.values() for lab in labs})
    print(f"Found {len(populations)} population(s) across {len(all_labels)} label(s): {all_labels}")
    if ignored:
        print(f"  ignored {len(ignored)} file(s) not matching the naming pattern")

    # every label present in the barcode files must have a BAM, or the extracted
    # population would be missing cells without saying so
    missing = [lab for lab in all_labels if lab not in bam_map]
    if missing:
        _fail(f"No --bam given for label(s) {missing}, which the barcode files require. "
              f"Every label needs one, otherwise the extracted populations would silently "
              f"be missing those cells. Provided: {sorted(bam_map)}")

    combo_ids = sorted(populations)
    if args.combo_ids:
        unknown = sorted(set(args.combo_ids) - set(combo_ids))
        if unknown:
            _fail(f"--combo-ids not found in {barcode_dir}: {unknown}")
        combo_ids = [c for c in combo_ids if c in set(args.combo_ids)]

    out_dir = os.path.abspath(args.output)
    script_dir = os.path.join(out_dir, "scripts")
    os.makedirs(script_dir, exist_ok=True)

    written, n_passes = [], 0
    for combo_id in combo_ids:
        per_label = [(lab, bam_map[lab], bc)
                     for lab, bc in sorted(populations[combo_id].items())]
        script_path = os.path.join(script_dir, f"combo_{combo_id}.sh")
        write_combo_script(script_path, combo_id, per_label, out_dir,
                           args.nproc, args.output_fragment)
        written.append(script_path)
        n_passes += len(per_label)

    runner = os.path.join(out_dir, "run_all.sh")
    write_runner(runner, written)

    print(f"\nWrote {len(written)} extraction scripts to: {script_dir}/")
    print(f"  {n_passes} sinto passes in total, one per label per population")
    print("\nNothing has been run. To extract one population:")
    print(f"    bash {written[0]}")
    print(f"To extract all {len(written)}:")
    print(f"    bash {runner}")


if __name__ == '__main__':
    main(argparse.Namespace())
