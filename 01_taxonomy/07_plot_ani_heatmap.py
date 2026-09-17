#!/usr/bin/env python3
"""
Turn fastANI's raw output (05_download_refs_and_fastani.slurm's
fastani_results.tsv) into a heatmap image.

fastANI's output has no header and is one row per query-reference pair
that cleared its internal detection floor (~ANI 75-80%; pairs too
distant to align at all are simply absent from the file, not zero):
    query_path  reference_path  ANI  matched_fragments  total_fragments

Usage (run this on a compute node, or the login node is fine too - this
is lightweight, not a real compute job):
    module load miniconda
    source activate /project/YOUR_PROJECT_HERE/envs/fastani_env
    python 07_plot_ani_heatmap.py \
        --input  fastani_out/fastani_results.tsv \
        --output fastani_out/ani_heatmap.png \
        --labels labels.tsv        # optional, see below

Optional --labels file: two columns (tab-separated, no header) mapping
an identifier substring to a display name, e.g.:
    Sample2               C. takakiae (Sample2)
    GCF_900129385.1       C. takakiae (ref)
    Sample4               Telluria sp. (Sample4)
    GCA_003484545.1       Telluria sp003484545 (ref)
Any query/reference whose file path contains that substring gets
relabeled; anything not matched falls back to the bare filename
(no directory, no extension).
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")  # no display on a compute/login node
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

# Sequential "blue" ramp, light -> dark, from the validated reference
# palette (dataviz skill, references/palette.md). One hue for one
# magnitude (ANI %) - no rainbow.
SEQUENTIAL_BLUE = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
    "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
    "#0d366b",
]
MISSING_COLOR = "#c3c2b7"   # muted gray - "not detected", not "zero"
GRIDLINE = "#e1e0d9"
INK = "#0b0b0b"
MUTED_INK = "#52514e"


def load_labels(path):
    mapping = []
    if not path:
        return mapping
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            mapping.append((parts[0], parts[1]))
    return mapping


def relabel(path, mapping):
    for needle, label in mapping:
        if needle in path:
            return label
    # fall back: bare filename, no directory, no extension
    base = os.path.basename(path)
    for ext in (".fna", ".fasta", ".fa", ".gz"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="fastANI output TSV (no header)")
    ap.add_argument("--output", required=True, help="output image path, e.g. ani_heatmap.png")
    ap.add_argument("--labels", default=None, help="optional TSV: substring<TAB>display_label")
    ap.add_argument("--annot", action="store_true", default=True, help="write ANI values in cells (default on)")
    ap.add_argument("--no-annot", dest="annot", action="store_false")
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()

    cols = ["query", "reference", "ani", "matched_frags", "total_frags"]
    df = pd.read_csv(args.input, sep="\t", header=None, names=cols)
    if df.empty:
        sys.exit(f"No rows read from {args.input} - is the path right and non-empty?")

    mapping = load_labels(args.labels)
    df["query_label"] = df["query"].apply(lambda p: relabel(p, mapping))
    df["ref_label"] = df["reference"].apply(lambda p: relabel(p, mapping))

    # Preserve first-seen order (matches your query_list.txt / ref_list.txt
    # order) rather than alphabetizing, so related references stay grouped.
    query_order = list(dict.fromkeys(df["query_label"]))
    ref_order = list(dict.fromkeys(df["ref_label"]))

    matrix = df.pivot_table(index="query_label", columns="ref_label", values="ani", aggfunc="mean")
    matrix = matrix.reindex(index=query_order, columns=ref_order)

    n_rows, n_cols = matrix.shape
    print(f"{n_rows} queries x {n_cols} references, {df.shape[0]} pairs with an ANI value "
          f"({n_rows * n_cols - df.shape[0]} pairs absent - below fastANI's detection floor).")

    cmap = LinearSegmentedColormap.from_list("ani_blue", SEQUENTIAL_BLUE, N=256)
    cmap.set_bad(MISSING_COLOR)

    fig_w = max(6.0, 0.55 * n_cols + 2.5)
    fig_h = max(3.5, 0.45 * n_rows + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("#fcfcfb")
    ax.set_facecolor("#fcfcfb")

    data = np.ma.masked_invalid(matrix.values)
    # Fix the scale to a sensible ANI range rather than autoscaling to
    # whatever happens to be in this run's data, so color is comparable
    # across re-runs with different reference sets.
    vmin, vmax = 75, 100
    im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(matrix.columns, rotation=90, ha="center", fontsize=8, color=MUTED_INK)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(matrix.index, fontsize=9, color=MUTED_INK)

    # Light gridlines between cells, not through them.
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color=GRIDLINE, linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if args.annot:
        for i in range(n_rows):
            for j in range(n_cols):
                v = matrix.values[i, j]
                if np.isnan(v):
                    continue
                # Flip label ink for readability on dark cells.
                frac = (v - vmin) / (vmax - vmin)
                txt_color = "#ffffff" if frac > 0.55 else INK
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6.5, color=txt_color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("ANI (%)", color=MUTED_INK)
    cbar.ax.tick_params(colors=MUTED_INK)

    ax.set_title("Average nucleotide identity (fastANI)", color=INK, fontsize=12, pad=12)
    ax.set_xlabel("Reference genome", color=MUTED_INK, fontsize=9)
    ax.set_ylabel("Query genome", color=MUTED_INK, fontsize=9)

    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi, facecolor=fig.get_facecolor())
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
