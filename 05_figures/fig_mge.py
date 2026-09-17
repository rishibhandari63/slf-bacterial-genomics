#!/usr/bin/env python3
"""Figures 8 and 9 - mobile genetic elements, resistance genes, defense systems."""
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

S = ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample8"]
LAB = ["S1  $\\it{A.\\ schindleri}$", "S2  $\\it{C.\\ takakiae}$",
       "S3  $\\it{S.\\ yabuuchiae}$", "S4  $\\it{Massilia}$ sp. nov.",
       "S5  $\\it{G.\\ mishrai}$", "S8  $\\it{Agromyces}$ sp. nov."]
NOVEL = [False, False, False, True, False, True]
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
SEQ = LinearSegmentedColormap.from_list(
    "seqblue", ["#f4f8fd", "#d3e4f7", "#a8c9ee", "#75a8e2", "#2a78d6", "#1b4f8f"])
INK, INK2 = "#0b0b0b", "#52514e"

sm = {r["Sample"]: r for r in csv.DictReader(open("mge_summary.csv"))}
g = lambda k: np.array([int(sm[s][k]) for s in S])

PRO, PLA = g("genomad_viruses"), g("genomad_plasmids")
AMR, STR = g("amr_AMR"), g("amr_STRESS")
DEF, DTY = g("defense_systems"), g("defense_types")

y = np.arange(len(S))


def styleax(ax, title, xlabel, labels=True):
    ax.set_yticks(y)
    ax.set_yticklabels(LAB if labels else [""] * len(S), fontsize=7.2)
    if labels:
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


# ============================ Figure 8 ============================
fig, (a, b, c) = plt.subplots(1, 3, figsize=(8.0, 2.9),
                              gridspec_kw={"width_ratios": [1.5, 1, 1]})
h = 0.36
a.barh(y - h / 2, PRO, h, color="#4a3aa7", linewidth=0, label="Prophages")
a.barh(y + h / 2, PLA, h, color="#1baf7a", linewidth=0, label="Plasmids")
for i in range(len(S)):
    if PRO[i]:
        a.text(PRO[i] + 0.1, i - h / 2, PRO[i], va="center", fontsize=6.8, color=INK2)
    if PLA[i]:
        a.text(PLA[i] + 0.1, i + h / 2, PLA[i], va="center", fontsize=6.8, color=INK2)
styleax(a, "a  Mobile genetic elements", "Elements detected")
a.set_xlim(0, max(PRO.max(), PLA.max()) * 1.25)
a.legend(fontsize=6.8, frameon=False, loc="lower right", handlelength=1.1)

b.barh(y - h / 2, AMR, h, color="#e34948", linewidth=0, label="AMR")
b.barh(y + h / 2, STR, h, color="#eda100", linewidth=0, label="Stress/metal")
for i in range(len(S)):
    if AMR[i]:
        b.text(AMR[i] + 0.12, i - h / 2, AMR[i], va="center", fontsize=6.8, color=INK2)
    if STR[i]:
        b.text(STR[i] + 0.12, i + h / 2, STR[i], va="center", fontsize=6.8, color=INK2)
styleax(b, "b  Resistance genes", "Genes detected", labels=False)
b.set_xlim(0, max(AMR.max(), STR.max()) * 1.3)
b.legend(fontsize=6.8, frameon=False, loc="lower right", handlelength=1.1)

c.barh(y - h / 2, DEF, h, color="#2a78d6", linewidth=0, label="Systems")
c.barh(y + h / 2, DTY, h, color="#e87ba4", linewidth=0, label="Distinct types")
for i in range(len(S)):
    if DEF[i]:
        c.text(DEF[i] + 0.3, i - h / 2, DEF[i], va="center", fontsize=6.8, color=INK2)
    if DTY[i]:
        c.text(DTY[i] + 0.3, i + h / 2, DTY[i], va="center", fontsize=6.8, color=INK2)
styleax(c, "c  Anti-phage defense", "Count", labels=False)
c.set_xlim(0, DEF.max() * 1.25)
c.legend(fontsize=6.8, frameon=False, loc="lower right", handlelength=1.1)

fig.tight_layout(w_pad=1.2)
fig.savefig("Fig8_MGE_defense.png", bbox_inches="tight")
fig.savefig("Fig8_MGE_defense.pdf", bbox_inches="tight")
plt.close(fig)

# ============================ Figure 9 ============================
rows = list(csv.reader(open("defense_system_counts.csv")))
hdr, data = rows[0], [r for r in rows[1:] if r]
idx = {s: hdr.index(s) for s in S}
types = [r[0] for r in data]
D = np.array([[int(r[idx[s]]) for s in S] for r in data])
order = np.argsort(-D.sum(axis=1))
types = [types[i] for i in order]
D = D[order]

fig, ax = plt.subplots(figsize=(4.3, 0.235 * len(types) + 1.7))
im = ax.imshow(D, cmap=SEQ, aspect="auto", vmin=0)
ax.set_xticks(range(len(S)))
ax.set_xticklabels(LAB, fontsize=7, rotation=42, ha="right", rotation_mode="anchor")
for t, nv in zip(ax.get_xticklabels(), NOVEL):
    if nv:
        t.set_fontweight("bold")
ax.set_yticks(range(len(types)))
ax.set_yticklabels(types, fontsize=6.9)
for i in range(D.shape[0]):
    for j in range(D.shape[1]):
        if D[i, j]:
            ax.text(j, i, D[i, j], ha="center", va="center", fontsize=6.4,
                    color="white" if D[i, j] > D.max() * 0.55 else INK2)
ax.set_xticks(np.arange(-.5, len(S), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(types), 1), minor=True)
ax.grid(which="minor", color="white", linewidth=1.4)
ax.tick_params(which="minor", length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, ticks=[0, 1, 2])
cb.set_label("Systems", fontsize=7, color=INK2)
cb.ax.tick_params(labelsize=6.5, length=2)
cb.outline.set_visible(False)
ax.set_title("Anti-phage defense system types", fontsize=9.2, color=INK,
             fontweight="bold", loc="left", pad=8)
fig.tight_layout()
fig.savefig("Fig9_defense_types.png", bbox_inches="tight")
fig.savefig("Fig9_defense_types.pdf", bbox_inches="tight")
print("written Fig8, Fig9")
print("defense types shown:", len(types))
