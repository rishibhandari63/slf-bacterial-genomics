#!/usr/bin/env python3
"""Figure 1 - genome, taxonomic and functional overview of the six isolates.

Consolidates the per-genome summary tables (genome statistics, dDDH to the
closest type strain, and functional/biosynthetic density) into one panel.
Colour identifies the genome and is held constant across all four panels, so
a reader tracks one isolate through every measure.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.linewidth": 0.6, "axes.edgecolor": "#8a8a85",
    "xtick.color": "#52514e", "ytick.color": "#52514e",
    "savefig.dpi": 400, "figure.dpi": 130,
})

S = ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample8"]
LAB = ["S1  $\\it{A.\\ schindleri}$", "S2  $\\it{C.\\ takakiae}$",
       "S3  $\\it{S.\\ yabuuchiae}$", "S4  $\\it{Massilia}$ sp. nov.",
       "S5  $\\it{G.\\ mishrai}$", "S8  $\\it{Agromyces}$ sp. nov."]
NOVEL = [False, False, False, True, False, True]
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
INK, INK2 = "#0b0b0b", "#52514e"

MB    = np.array([3.281, 4.825, 4.326, 5.706, 3.540, 4.281])
GC    = np.array([42.5, 37.1, 66.0, 65.9, 59.5, 71.3])
DDDH  = np.array([73.0, 90.9, 73.6, 51.3, 88.4, 33.2])
LO    = np.array([70.0, 88.7, 70.6, 48.6, 85.9, 30.8])
HI    = np.array([75.8, 92.7, 76.4, 53.9, 90.5, 35.7])
CAZMB = np.array([15.5, 34.6, 42.3, 24.4, 26.0, 43.7])
BGCMB = np.array([1.22, 0.83, 0.92, 1.58, 1.98, 1.40])

y = np.arange(len(S))


def base(ax, title, xlabel, show_labels):
    ax.set_yticks(y)
    ax.set_yticklabels(LAB if show_labels else [""] * len(S), fontsize=7.4)
    if show_labels:
        for t, nv in zip(ax.get_yticklabels(), NOVEL):
            if nv:
                t.set_fontweight("bold")
    ax.invert_yaxis()
    ax.set_xlabel(xlabel, fontsize=7.8, color=INK2)
    ax.set_title(title, fontsize=8.8, color=INK, fontweight="bold", loc="left", pad=6)
    ax.grid(axis="x", color="#e8e8e4", linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(length=2.5, labelsize=7.2)


fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.4))
(a, b), (c, d) = axes

# ---- a: genome size, with G+C annotated ----
a.barh(y, MB, 0.62, color=CAT, linewidth=0)
for i, (v, g) in enumerate(zip(MB, GC)):
    a.text(v + 0.09, i, f"{v:.2f} Mb   {g:.1f}% G+C", va="center",
           fontsize=6.8, color=INK2)
base(a, "a  Assembly size and G+C content", "Genome size (Mb)", True)
a.set_xlim(0, MB.max() * 1.62)

# ---- b: dDDH with confidence intervals and the 70% species threshold ----
err = np.vstack([DDDH - LO, HI - DDDH])
b.barh(y, DDDH, 0.62, color=CAT, linewidth=0,
       xerr=err, error_kw=dict(ecolor="#52514e", elinewidth=0.9, capsize=2.2))
b.axvline(70, color="#e34948", linestyle="--", linewidth=1.1, zorder=3)
b.text(70, -0.92, "70% species\nthreshold", fontsize=6.5, color="#e34948",
       ha="center", va="bottom", linespacing=1.15)
for i, v in enumerate(DDDH):
    b.text(min(v + (HI - DDDH)[i] + 2.0, 101), i, f"{v:.1f}", va="center",
           fontsize=6.8, color=INK2)
base(b, "b  dDDH to closest type strain", "dDDH, formula d4 (%)", False)
b.set_xlim(0, 112)
b.set_ylim(len(S) - 0.4, -1.35)

# ---- c: CAZyme density ----
c.barh(y, CAZMB, 0.62, color=CAT, linewidth=0)
for i, v in enumerate(CAZMB):
    c.text(v + 0.7, i, f"{v:.1f}", va="center", fontsize=6.8, color=INK2)
base(c, "c  Carbohydrate-active enzymes", "CAZymes per Mb", True)
c.set_xlim(0, CAZMB.max() * 1.16)

# ---- d: BGC density ----
d.barh(y, BGCMB, 0.62, color=CAT, linewidth=0)
for i, v in enumerate(BGCMB):
    d.text(v + 0.035, i, f"{v:.2f}", va="center", fontsize=6.8, color=INK2)
base(d, "d  Biosynthetic gene clusters", "BGC regions per Mb", False)
d.set_xlim(0, BGCMB.max() * 1.18)

fig.tight_layout(w_pad=1.4, h_pad=2.0)
fig.savefig("Fig1_genome_overview.png", bbox_inches="tight")
fig.savefig("Fig1_genome_overview.pdf", bbox_inches="tight")
print("written Fig1_genome_overview")
