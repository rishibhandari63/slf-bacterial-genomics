#!/usr/bin/env python3
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.linewidth": 0.6, "axes.edgecolor": "#8a8a85",
    "xtick.color": "#52514e", "ytick.color": "#52514e",
    "savefig.dpi": 400, "figure.dpi": 130,
})

SAMPLES = ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample8"]
LAB = {"Sample1": "S1  A. schindleri", "Sample2": "S2  C. takakiae",
       "Sample3": "S3  S. yabuuchiae", "Sample4": "S4  Massilia sp. nov.",
       "Sample5": "S5  G. mishrai", "Sample8": "S8  Agromyces sp. nov."}
MB = {"Sample1": 3.281, "Sample2": 4.825, "Sample3": 4.326,
      "Sample4": 5.706, "Sample5": 3.540, "Sample8": 4.281}
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
SEQ = LinearSegmentedColormap.from_list(
    "seqblue", ["#f4f8fd", "#d3e4f7", "#a8c9ee", "#75a8e2", "#2a78d6", "#1b4f8f"])
INK, INK2 = "#0b0b0b", "#52514e"

rows = list(csv.reader(open("bgc_type_counts.csv")))
hdr, data = rows[0], [r for r in rows[1:] if r]
idx = {s: hdr.index(s) for s in SAMPLES}
types = [r[0] for r in data]
T = np.array([[int(r[idx[s]]) for s in SAMPLES] for r in data])
order = np.argsort(-T.sum(axis=1))
types = [types[i] for i in order]
T = T[order]

# IMPORTANT: do NOT sum the type matrix to get region counts. Hybrid regions
# carry several product types and are counted once per type, so the column
# sums exceed the true number of regions (e.g. Sample8: 8 type-assignments
# across 6 regions). Region counts come from bgc_summary.csv.
srows = list(csv.reader(open("bgc_summary.csv")))
shdr = srows[0]
sreg = {r[0]: int(r[shdr.index("Total BGC regions")]) for r in srows[1:] if r}
totals = np.array([sreg[s] for s in SAMPLES])
type_assignments = T.sum(axis=0)
dens = totals / np.array([MB[s] for s in SAMPLES])

fig, (ax1, ax2) = plt.subplots(
    1, 2, figsize=(7.6, 4.5), gridspec_kw={"width_ratios": [1.75, 1]})

im = ax1.imshow(T, cmap=SEQ, aspect="auto", vmin=0)
ax1.set_xticks(range(len(SAMPLES)))
ax1.set_xticklabels([LAB[s] for s in SAMPLES], fontsize=7, rotation=40,
                    ha="right", rotation_mode="anchor")
for t in ax1.get_xticklabels():
    if "nov." in t.get_text():
        t.set_fontweight("bold")
ax1.set_yticks(range(len(types)))
ax1.set_yticklabels(types, fontsize=7.2)
for i in range(T.shape[0]):
    for j in range(T.shape[1]):
        if T[i, j]:
            ax1.text(j, i, T[i, j], ha="center", va="center", fontsize=6.8,
                     color="white" if T[i, j] > T.max() * 0.55 else INK2)
ax1.set_xticks(np.arange(-.5, len(SAMPLES), 1), minor=True)
ax1.set_yticks(np.arange(-.5, len(types), 1), minor=True)
ax1.grid(which="minor", color="white", linewidth=1.5)
ax1.tick_params(which="minor", length=0)
for s in ax1.spines.values():
    s.set_visible(False)
ax1.set_title("a  BGC product types", fontsize=9, color=INK,
              fontweight="bold", loc="left", pad=8)
cb = fig.colorbar(im, ax=ax1, fraction=0.03, pad=0.02)
cb.set_label("Regions", fontsize=7, color=INK2)
cb.ax.tick_params(labelsize=6.5, length=2)
cb.outline.set_visible(False)

y = np.arange(len(SAMPLES))
ax2.barh(y, dens, 0.6, color=CAT, linewidth=0)
for i, (d, tot) in enumerate(zip(dens, totals)):
    ax2.text(d + 0.05, i, f"{d:.2f}  (n={tot})", va="center", fontsize=6.9, color=INK2)
ax2.set_yticks(y)
ax2.set_yticklabels([LAB[s] for s in SAMPLES], fontsize=7)
for t in ax2.get_yticklabels():
    if "nov." in t.get_text():
        t.set_fontweight("bold")
ax2.invert_yaxis()
ax2.set_xlabel("BGC regions per Mb", fontsize=8, color=INK2)
ax2.set_xlim(0, dens.max() * 1.42)
ax2.grid(axis="x", color="#e8e8e4", linewidth=0.6)
ax2.set_axisbelow(True)
for s in ("top", "right"):
    ax2.spines[s].set_visible(False)
ax2.tick_params(length=2.5)
ax2.set_title("b  BGC density", fontsize=9, color=INK, fontweight="bold",
              loc="left", pad=8)

fig.tight_layout()
fig.savefig("Fig7_BGC.png", bbox_inches="tight")
fig.savefig("Fig7_BGC.pdf", bbox_inches="tight")
print("regions:", dict(zip(SAMPLES, totals.tolist())))
print("type-assignments (hybrids counted per type):",
      dict(zip(SAMPLES, type_assignments.tolist())))
print("per Mb:", {s: round(float(d), 2) for s, d in zip(SAMPLES, dens)})
