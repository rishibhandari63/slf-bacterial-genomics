#!/usr/bin/env python3
"""Composite Figure 1 (7.5 x 8 in): (a) genome size and G+C, (b) bac120 phylogeny,
(c) dDDH to closest type strain, (d) ANIb alignment coverage matrix.

Panels a and c are drawn here. Panels b and d are embedded as VECTOR from the
supplied PDFs with show_pdf_page(clip=...), with all labelling redrawn at
figure-appropriate sizes because the source labels scale to ~2 pt.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
import pymupdf

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7,
    "axes.linewidth": 0.5, "axes.edgecolor": "#8a8a85",
    "xtick.color": "#52514e", "ytick.color": "#52514e",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

W_PT, H_PT = 540.0, 576.0          # 7.5 x 8 inch
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e8e8e4"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]

# ---------------------------------------------------------------- panel data
LAB_TAG = ["S1", "S2", "S3", "S4", "S5", "S8"]
LAB_SP  = ["A. schindleri", "C. takakiae", "S. yabuuchiae", "Massilia", "G. mishrai", "Agromyces"]
LAB_SUF = ["", "", "", " sp. nov.", "", " sp. nov."]
NOVEL = [False, False, False, True, False, True]
MB   = np.array([3.281, 4.825, 4.326, 5.706, 3.540, 4.281])
GC   = np.array([42.5, 37.1, 66.0, 65.9, 59.5, 71.3])
DDDH = np.array([73.0, 90.9, 73.6, 51.3, 88.4, 33.2])
LO   = np.array([70.0, 88.7, 70.6, 48.6, 85.9, 30.8])
HI   = np.array([75.8, 92.7, 76.4, 53.9, 90.5, 35.7])
y = np.arange(6)

# ------------------------------------------------------------ panel rects (pt, top-left origin)
A_R = (6, 6, 262, 206)
C_R = (288, 6, 536, 206)
B_R = (6, 222, 262, 574)
D_R = (288, 222, 536, 574)


def ax_from_rect(fig, r, pad=(0, 0, 0, 0)):
    """rect in pt (x0,y0top,x1,y1top) -> matplotlib axes in figure fractions."""
    x0, y0, x1, y1 = r
    x0 += pad[0]; y0 += pad[1]; x1 -= pad[2]; y1 -= pad[3]
    return fig.add_axes([x0 / W_PT, 1 - y1 / H_PT, (x1 - x0) / W_PT, (y1 - y0) / H_PT])


def panel_letter(fig, r, letter):
    fig.text(r[0] / W_PT, 1 - (r[1] + 7) / H_PT, letter, fontsize=9.5,
             fontweight="bold", color=INK, ha="left", va="center")


def tidy(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="x", color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(length=2, labelsize=6.4, pad=1.5)
    ax.invert_yaxis()


fig = plt.figure(figsize=(W_PT / 72, H_PT / 72))
fig.patch.set_alpha(0.0)

# ----------------------------------------------------------------- panel (a)
a = ax_from_rect(fig, A_R, pad=(94, 14, 8, 24))
a.barh(y, MB, color=CAT, height=0.66, linewidth=0)
for i, (m, g) in enumerate(zip(MB, GC)):
    a.text(m + 0.09, i, f"{g:.1f}% G+C", va="center", fontsize=5.9, color=INK2)
a.set_yticks(y)
a.set_yticklabels([])
a.tick_params(axis="y", length=0)
PENDING = []
for i in range(6):
    fw = "bold" if NOVEL[i] else "normal"
    a.annotate(LAB_TAG[i], xy=(0, i), xycoords=("axes fraction", "data"),
               xytext=(-88, 0), textcoords="offset points", fontsize=6.2,
               fontweight=fw, color=INK, va="center", ha="left", annotation_clip=False)
    sp_txt = a.annotate(LAB_SP[i], xy=(0, i), xycoords=("axes fraction", "data"),
               xytext=(-74, 0), textcoords="offset points", fontsize=6.2, style="italic",
               fontweight=fw, color=INK, va="center", ha="left", annotation_clip=False)
    if LAB_SUF[i]:
        PENDING.append((i, fw, sp_txt))
a.set_xlim(0, 7.4)
a.set_xlabel("Assembly size (Mb)", fontsize=6.8, color=INK2, labelpad=1.5)
tidy(a)
a.tick_params(axis="y", length=0)
fig.set_dpi(72)
_r = fig.canvas.get_renderer()
for i, fw, t in PENDING:
    wpt = t.get_window_extent(renderer=_r).width
    a.annotate("sp. nov.", xy=(0, i), xycoords=("axes fraction", "data"),
               xytext=(-74 + wpt + 2.6, 0), textcoords="offset points", fontsize=6.2,
               fontweight=fw, color=INK, va="center", ha="left", annotation_clip=False)
panel_letter(fig, A_R, "a")

# ----------------------------------------------------------------- panel (c)
c = ax_from_rect(fig, C_R, pad=(16, 14, 8, 24))
c.barh(y, DDDH, color=CAT, height=0.66, linewidth=0)
c.errorbar(DDDH, y, xerr=[DDDH - LO, HI - DDDH], fmt="none",
           ecolor=INK2, elinewidth=0.6, capsize=1.6, capthick=0.6)
c.axvline(70, color="#e34948", linestyle="--", linewidth=0.8)
c.text(70.8, -0.72, "70% species threshold", fontsize=5.7, color="#e34948",
       ha="left", va="center")
c.set_yticks(y)
c.set_yticklabels([])
c.set_xlim(0, 100)
c.set_xlabel("dDDH to closest type strain (%)", fontsize=6.8, color=INK2, labelpad=1.5)
tidy(c)
c.tick_params(axis="y", length=0)
panel_letter(fig, C_R, "c")

# ----------------------------------------------------------------- panel (b) labels
TIPS = [
    ("outgroup", "C. trachomatis (outgroup)", None),
    ("Chryseobacterium", "C. sp. SP36", None),
    ("Chryseobacterium", "C. shandongense H5143", None),
    ("Chryseobacterium", "C. hispalense DSM 25574", None),
    ("Chryseobacterium", "C. nepalense AC3", None),
    ("Chryseobacterium", "C. profundimaris DSM 28214", None),
    ("Chryseobacterium", "C. takakiae DSM 26898", None),
    ("Chryseobacterium", "Sample2", 1),
    ("Glutamicibacter", "G. nicotianae NBRC 14234", None),
    ("Glutamicibacter", "Arthrobacter sp. NIO-1057", None),
    ("Glutamicibacter", "G. halophytocola KLBMP 5180", None),
    ("Glutamicibacter", "G. sp. JL-03c", None),
    ("Glutamicibacter", "G. mishrai S5-52", None),
    ("Glutamicibacter", "Sample5", 4),
    ("Agromyces", "A. salentinus JCM 14323", None),
    ("Agromyces", "A. terreus JCM 14581", None),
    ("Agromyces", "A. sp. CF514", None),
    ("Agromyces", "Sample8", 5),
    ("Agromyces", "A. sp. Leaf222", None),
    ("Agromyces", "A. aureus AR33", None),
    ("Agromyces", "A. allii JCM 13584", None),
    ("Sphingomonas", "S. sanguinis NBRC 13937", None),
    ("Sphingomonas", "S. sanguinis NP2-R2", None),
    ("Sphingomonas", "S. excrementigallinarum 1562", None),
    ("Sphingomonas", "S. yabuuchiae DSM 14562", None),
    ("Sphingomonas", "Sample3", 2),
    ("Sphingomonas", "S. parapaucimobilis NBRC 15100", None),
    ("Sphingomonas", "S. yabuuchiae NS355", None),
    ("Acinetobacter", "A. variabilis NIPH 2171", None),
    ("Acinetobacter", "A. lwoffii H7", None),
    ("Acinetobacter", "A. idrijaensis MII", None),
    ("Acinetobacter", "A. sp. LoGeW2-3", None),
    ("Acinetobacter", "A. sp. ASP199", None),
    ("Acinetobacter", "A. schindleri CIP 107287", None),
    ("Acinetobacter", "Sample1", 0),
    ("Massilia", "M. sp. YIMB02769", None),
    ("Massilia", "M. brevitalea CGMCC 1.10731", None),
    ("Massilia", "M. sp. IC2-278", None),
    ("Massilia", "M. sp. UBA11196", None),
    ("Massilia", "M. sp. UBA2709", None),
    ("Massilia", "M. sp. UBA9086", None),
    ("Massilia", "Sample4", 3),
]
SAMPLE_NAME = {"Sample1": "Sample1  (A. schindleri)", "Sample2": "Sample2  (C. takakiae)",
               "Sample3": "Sample3  (S. yabuuchiae)", "Sample4": "Sample4  sp. nov.",
               "Sample5": "Sample5  (G. mishrai)", "Sample8": "Sample8  sp. nov."}

DEND_X0, DEND_W = 17.0, 130.0           # dendrogram band, absolute page pt
LAB_X = DEND_X0 + DEND_W + 4.0
BR_X = 239.0                             # genus bracket x
TREE_TOP, TREE_BOT = 244.0, 572.0        # tip band (pt, top-left origin)

b = ax_from_rect(fig, B_R)
b.set_xlim(B_R[0], B_R[2]); b.set_ylim(B_R[3], B_R[1])
b.axis("off")
pitch = (TREE_BOT - TREE_TOP) / 42.0
tip_y = [TREE_TOP + (i + 0.5) * pitch for i in range(42)]

for (genus, txt, si), yy in zip(TIPS, tip_y):
    if txt.startswith("Sample"):
        b.text(LAB_X, yy, SAMPLE_NAME[txt], fontsize=5.0, style="normal",
               fontweight="bold", color=CAT[si], va="center", ha="left")
    elif genus == "outgroup":
        b.text(LAB_X, yy, txt, fontsize=5.0, style="italic", color=INK2,
               va="center", ha="left")
    else:
        b.text(LAB_X, yy, txt, fontsize=5.0, style="italic", color="#3a3a37",
               va="center", ha="left")

# genus brackets
blocks, start = [], 1
for i in range(1, 43):
    if i == 42 or TIPS[i][0] != TIPS[start][0]:
        blocks.append((start, i - 1, TIPS[start][0])); start = i
for i0, i1, g in blocks:
    y0, y1 = tip_y[i0] - pitch * 0.33, tip_y[i1] + pitch * 0.33
    b.plot([BR_X, BR_X], [y0, y1], color="#9a9a95", linewidth=0.7, clip_on=False)
    b.plot([BR_X, BR_X - 2.2], [y0, y0], color="#9a9a95", linewidth=0.7, clip_on=False)
    b.plot([BR_X, BR_X - 2.2], [y1, y1], color="#9a9a95", linewidth=0.7, clip_on=False)
    b.text(BR_X + 2.4, (y0 + y1) / 2, g, fontsize=5.4, style="italic", color=INK2,
           rotation=90, va="center", ha="center", rotation_mode="anchor")

# scale bar: 58.88 pt = 0.1 subs/site in the source, source topology width 725.5 pt
SRC_W = 787.3 - 60.5
subs_per_pt = 0.1 / (58.88 * DEND_W / SRC_W)
bar_val = 0.2
bar_len = bar_val / subs_per_pt
b.plot([DEND_X0, DEND_X0 + bar_len], [TREE_TOP - 11, TREE_TOP - 11],
       color=INK2, linewidth=0.8, clip_on=False)
for xx in (DEND_X0, DEND_X0 + bar_len):
    b.plot([xx, xx], [TREE_TOP - 13.2, TREE_TOP - 8.8], color=INK2, linewidth=0.8, clip_on=False)
b.text(DEND_X0 + bar_len + 3, TREE_TOP - 11, f"{bar_val} subs/site", fontsize=5.4,
       color=INK2, va="center", ha="left")
panel_letter(fig, B_R, "b")

# ----------------------------------------------------------------- panel (d) annotation
ORDER_BLOCKS = [(0, 6, "Chryseobacterium"), (7, 13, "Acinetobacter"),
                (14, 19, "Glutamicibacter"), (20, 26, "Agromyces"),
                (27, 33, "Sphingomonas"), (34, 40, "Massilia")]
SAMPLE_AT = {2: "S2", 10: "S1", 19: "S5", 22: "S8", 33: "S3", 38: "S4"}

MAT_X0, MAT_X1 = 306.0, 536.0
MAT_Y0, MAT_Y1 = 262.0, 492.0
d = ax_from_rect(fig, D_R)
d.set_xlim(D_R[0], D_R[2]); d.set_ylim(D_R[3], D_R[1])
d.axis("off")
n = 41
cw = (MAT_X1 - MAT_X0) / n

d.add_patch(plt.Rectangle((MAT_X0, MAT_Y0), MAT_X1 - MAT_X0, MAT_Y1 - MAT_Y0,
                          fill=False, edgecolor="#8a8a85", linewidth=0.5, zorder=5))
for i0, i1, g in ORDER_BLOCKS:
    a0, a1 = MAT_X0 + i0 * cw, MAT_X0 + (i1 + 1) * cw
    b0, b1 = MAT_Y0 + i0 * cw, MAT_Y0 + (i1 + 1) * cw
    d.add_patch(plt.Rectangle((a0, b0), a1 - a0, b1 - b0, fill=False,
                              edgecolor="#ffffff", linewidth=0.9, zorder=6))
    d.text(MAT_X0 - 14.0, (b0 + b1) / 2, g, fontsize=4.4, style="italic", color=INK2,
           rotation=90, va="center", ha="center", rotation_mode="anchor", zorder=7)
TAGS = "S1 S2 S3 S4 S5 S8".split()
for idx, tag in SAMPLE_AT.items():
    col = CAT[TAGS.index(tag)]
    yy = MAT_Y0 + (idx + 0.5) * cw
    d.plot([MAT_X0 - 2.2], [yy], marker=">", markersize=2.4, color=col,
           clip_on=False, zorder=8, linestyle="none")
    d.text(MAT_X0 - 5.0, yy, tag, fontsize=4.6, fontweight="bold", color=col,
           va="center", ha="right", zorder=8)

# colour bar reproducing pyani's blue-white-red 0-1 scale
cmap = LinearSegmentedColormap.from_list("pyani_bwr", ["#0000ff", "#ffffff", "#ff0000"])
cb_ax = fig.add_axes([(MAT_X0 + 45) / W_PT, 1 - 528.0 / H_PT, 140 / W_PT, 7.0 / H_PT])
cb = fig.colorbar(ScalarMappable(norm=Normalize(0, 1), cmap=cmap), cax=cb_ax,
                  orientation="horizontal")
cb.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
cb.ax.tick_params(labelsize=5.4, length=1.8, pad=1.2, color=INK2)
cb.outline.set_linewidth(0.4); cb.outline.set_edgecolor("#8a8a85")
cb.set_label("ANIb alignment coverage (fraction of genome aligned)",
             fontsize=5.8, color=INK2, labelpad=1.8)
panel_letter(fig, D_R, "d")

fig.savefig("overlay.pdf", transparent=True, dpi=600)
plt.close(fig)

# ------------------------------------------------------------------ composite
out = pymupdf.open()
page = out.new_page(width=W_PT, height=H_PT)
page.draw_rect(pymupdf.Rect(0, 0, W_PT, H_PT), color=None, fill=(1, 1, 1))

tree = pymupdf.open("tree_named.pdf")
TIP0, TIPN, PITCH_SRC = 53.0, 940.9, 21.66
CY0 = 49.6                                   # below the source "Tree scale" text
CY1 = TIPN + PITCH_SRC / 2                   # half a pitch past the last tip
pitch_t = (TREE_BOT - TREE_TOP) / 42.0
th = (CY1 - CY0) * pitch_t / PITCH_SRC       # so one source pitch == one target pitch
ty0 = TREE_BOT - th
clip_tree = pymupdf.Rect(60.5, CY0, 787.3, CY1)
page.show_pdf_page(pymupdf.Rect(DEND_X0, ty0, DEND_X0 + DEND_W, TREE_BOT),
                   tree, 0, clip=clip_tree, keep_proportion=False)

# strip the source text: the per-cell numbers are illegible at panel size and
# carry a Type 3 font that some journals reject
_c = pymupdf.open("cov.pdf"); _cp = _c[0]
_cp.add_redact_annot(_cp.rect)
try:
    _cp.apply_redactions(images=pymupdf.PDF_REDACT_IMAGE_NONE,
                         graphics=pymupdf.PDF_REDACT_LINE_ART_NONE,
                         text=pymupdf.PDF_REDACT_TEXT_REMOVE)
except Exception:
    _cp.apply_redactions()
_c.save("cov_notext.pdf", garbage=4, deflate=True)
cov = pymupdf.open("cov_notext.pdf")
page.show_pdf_page(pymupdf.Rect(MAT_X0, MAT_Y0, MAT_X1, MAT_Y1), cov, 0,
                   clip=pymupdf.Rect(724.9, 726.1, 3131.0, 3132.0), keep_proportion=False)

ov = pymupdf.open("overlay.pdf")
page.show_pdf_page(pymupdf.Rect(0, 0, W_PT, H_PT), ov, 0, keep_proportion=False)

out.save("Figure_1_composite.pdf", garbage=4, deflate=True)
print("written Figure_1_composite.pdf  %.2f x %.2f in" % (W_PT / 72, H_PT / 72))
