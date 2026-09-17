#!/usr/bin/env python3
"""Figure 10 - genus-specific gene content of the two candidate novel species."""
import csv, re
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

INFO = {
    "Sample8": dict(prefix="LGGIJPPN", n_prot=3841, n_ref=152,
                    lab="S8  $\\it{Agromyces}$ sp. nov.", col="#008300"),
    "Sample4": dict(prefix="GJCANLFC", n_prot=4919, n_ref=30,
                    lab="S4  $\\it{Massilia}$ sp. nov.", col="#eda100"),
}
ORDER = ["Sample8", "Sample4"]
INK, INK2 = "#0b0b0b", "#52514e"
C_ABS, C_DIV = "#4a3aa7", "#e87ba4"     # absent / divergent


def load(sample):
    d = INFO[sample]
    ident = []
    for r in csv.DictReader(open(f"{sample}_best_hit_identity.csv")):
        if r["Has homolog"] == "yes" and r["Best hit % identity"]:
            ident.append(float(r["Best hit % identity"]))
    sets = {}
    for tag in ("absent_from_genus", "divergent"):
        idx = []
        for r in csv.DictReader(open(f"{sample}_{tag}.csv")):
            m = re.match(rf"{d['prefix']}_(\d+)", r["Protein"])
            if m:
                idx.append(int(m.group(1)))
        sets[tag] = sorted(idx)
    cog = {}
    for r in csv.DictReader(open(f"{sample}_divergent_cog.csv")):
        if r["COG category"] != "-":
            cog[r["COG category"]] = (r["Description"], int(r["Proteins"]))
    return ident, sets, cog


data = {s: load(s) for s in ORDER}

fig = plt.figure(figsize=(7.6, 6.6))
gs = fig.add_gridspec(3, 2, height_ratios=[1.05, 0.72, 1.15], hspace=0.62, wspace=0.22)

# ---------------- a: best-hit identity distributions ----------------
for k, s in enumerate(ORDER):
    ax = fig.add_subplot(gs[0, k])
    ident, sets, _ = data[s]
    d = INFO[s]
    ax.hist(ident, bins=np.arange(20, 101, 2), color=d["col"], linewidth=0)
    ax.axvline(40, color="#e34948", linestyle="--", linewidth=1.0)
    ax.text(41, ax.get_ylim()[1] * 0.93, "40%", fontsize=6.6, color="#e34948")
    med = float(np.median(ident))
    ax.axvline(med, color=INK2, linestyle=":", linewidth=1.0)
    ax.text(med - 1.5, ax.get_ylim()[1] * 0.72, f"median {med:.1f}%", fontsize=6.6,
            color=INK2, ha="right")
    ax.set_xlabel("Best-hit identity to any congener (%)", fontsize=7.6, color=INK2)
    if k == 0:
        ax.set_ylabel("Proteins", fontsize=7.8, color=INK2)
    ax.set_title(f"{'ab'[k]}  {d['lab']}", fontsize=8.6, color=INK,
                 fontweight="bold", loc="left", pad=6)
    ax.set_xlim(20, 100)
    ax.grid(axis="y", color="#e8e8e4", linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(length=2.5, labelsize=7)
    n_abs, n_div = len(sets["absent_from_genus"]), len(sets["divergent"])
    ax.text(0.04, 0.95, f"{n_abs} no homolog\n{n_div} <40% identity\nn={d['n_ref']} congener genomes",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.6, color=INK2,
            linespacing=1.55)

# ---------------- c: positional distribution along the genome ----------------
axm = fig.add_subplot(gs[1, :])
for k, s in enumerate(ORDER):
    _, sets, _ = data[s]
    d = INFO[s]
    y = 1 - k
    axm.hlines(y, 0, d["n_prot"], color="#e0e0dc", linewidth=5, zorder=1)
    for tag, colr, off in (("absent_from_genus", C_ABS, 0.16),
                           ("divergent", C_DIV, -0.16)):
        xs = sets[tag]
        axm.vlines(xs, y + (0.055 if off > 0 else -0.20), y + (0.20 if off > 0 else -0.055),
                   color=colr, linewidth=1.0, zorder=3)
    axm.text(-d["n_prot"] * 0.012, y, d["lab"], ha="right", va="center",
             fontsize=7.4, fontweight="bold")
axm.set_ylim(-0.65, 1.65)
axm.set_xlim(-950, 5050)
axm.set_xlabel("Gene index along the genome", fontsize=7.8, color=INK2)
axm.set_yticks([])
for sp in ("top", "right", "left"):
    axm.spines[sp].set_visible(False)
axm.tick_params(length=2.5, labelsize=7)
axm.set_title("c  Genomic distribution of genus-specific genes", fontsize=8.6,
              color=INK, fontweight="bold", loc="left", pad=6)
h = [plt.Line2D([0], [0], color=C_ABS, lw=2), plt.Line2D([0], [0], color=C_DIV, lw=2)]
axm.legend(h, ["no homolog in genus", "best hit <40% identity"], fontsize=6.8,
           frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.30),
           handlelength=1.2)

# ---------------- d: COG profile of the divergent proteins ----------------
axc = fig.add_subplot(gs[2, :])
cats = []
for s in ORDER:
    cats += list(data[s][2].keys())
cats = sorted(set(cats), key=lambda c: -sum(data[s][2].get(c, ("", 0))[1] for s in ORDER))
x = np.arange(len(cats))
w = 0.38
for k, s in enumerate(ORDER):
    cog = data[s][2]
    vals = [cog.get(c, ("", 0))[1] for c in cats]
    axc.bar(x + (k - 0.5) * w, vals, w * 0.9, color=INFO[s]["col"], linewidth=0,
            label=INFO[s]["lab"])
    for xi, v in zip(x + (k - 0.5) * w, vals):
        if v:
            axc.text(xi, v + 0.12, v, ha="center", fontsize=6.3, color=INK2)
axc.set_xticks(x)
axc.set_xticklabels(cats, fontsize=7.6)
axc.set_ylabel("Divergent proteins", fontsize=7.8, color=INK2)
axc.set_title("d  COG categories of highly divergent proteins", fontsize=8.6,
              color=INK, fontweight="bold", loc="left", pad=6)
axc.grid(axis="y", color="#e8e8e4", linewidth=0.6)
axc.set_axisbelow(True)
for sp in ("top", "right"):
    axc.spines[sp].set_visible(False)
axc.tick_params(length=2.5, labelsize=7)
axc.legend(fontsize=6.9, frameon=False, loc="upper right", handlelength=1.1)
key = {"L":"Replication/recombination","K":"Transcription","U":"Secretion/trafficking",
       "S":"Function unknown","M":"Cell wall/membrane","P":"Inorganic ion transport",
       "O":"Post-translational modification","N":"Cell motility","T":"Signal transduction",
       "C":"Energy production","V":"Defense","H":"Coenzyme transport","J":"Translation",
       "G":"Carbohydrate transport","D":"Cell cycle","E":"Amino acid transport",
       "I":"Lipid transport","Q":"Secondary metabolites"}
axc.text(0, -0.42, "  ".join(f"{c}, {key.get(c,'')}" for c in cats[:9]),
         transform=axc.transAxes, fontsize=6.0, color=INK2)
axc.text(0, -0.55, "  ".join(f"{c}, {key.get(c,'')}" for c in cats[9:]),
         transform=axc.transAxes, fontsize=6.0, color=INK2)

fig.savefig("Fig10_unique_genes.png", bbox_inches="tight")
fig.savefig("Fig10_unique_genes.pdf", bbox_inches="tight")
print("written Fig10")
for s in ORDER:
    ident, sets, cog = data[s]
    print(f"  {s}: {len(ident)} with homolog, {len(sets['absent_from_genus'])} absent, "
          f"{len(sets['divergent'])} divergent, {len(cog)} COG categories")
