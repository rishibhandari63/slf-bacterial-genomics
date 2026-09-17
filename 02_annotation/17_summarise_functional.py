#!/usr/bin/env python3
"""
Summarise eggNOG-mapper + dbCAN output into tidy CSV tables.

Deliberately uses ONLY the Python standard library - no pandas, no matplotlib -
so it runs with whatever python3 is on the Ceres login node without touching
any conda environment. Figures are built separately from these CSVs.

Usage (login node is fine, this is lightweight). Type it on ONE line -
a "\" continuation followed by a space instead of a newline silently
corrupts the paths:

    python3 17_summarise_functional.py <annotation_dir> <prokka_dir>

e.g.
    python3 17_summarise_functional.py /project/kanglab/rishi.bhandari/funct_annotation_out /project/kanglab/rishi.bhandari/PROKKA_DIR/prokka_out

Writes into <outdir>/summary/:
    annotation_summary.csv   - per genome: proteins, annotated, COG-assigned, KEGG KO, CAZy
    cog_counts.csv           - COG category x genome count matrix (+ category descriptions)
    cog_percent.csv          - same, as % of COG-assigned proteins in that genome
    kegg_pathway_counts.csv  - KEGG ko pathway x genome counts (top pathways)
    cazyme_class_counts.csv  - CAZyme class (GH/GT/PL/CE/AA/CBM) x genome
    cazyme_family_counts.csv - CAZyme family x genome
"""
import csv
import os
import re
import sys
from collections import Counter, defaultdict

COG_DESC = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion, vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "-": "Unassigned",
}
CAZY_CLASSES = ["GH", "GT", "PL", "CE", "AA", "CBM"]


def count_faa(path):
    """Count sequences in the Prokka .faa that produced this annotation, if findable."""
    if not path or not os.path.isfile(path):
        return None
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def parse_emapper(path):
    """Yield dict rows from an .emapper.annotations file."""
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##") or not line.strip():
                continue
            if line.startswith("#query"):
                header = line.lstrip("#").split("\t")
                header[0] = "query"
                continue
            if header is None:
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            yield dict(zip(header, parts))


def parse_dbcan(sample_dir):
    """
    Return Counter of CAZyme families for one sample.

    dbCAN's output filename has changed across versions, so look for any
    overview-style table rather than hardcoding one name. A protein is counted
    as a CAZyme when at least 2 of the 3 tools agree (the standard dbCAN
    convention); if no '#ofTools'-style column exists we fall back to counting
    any row with a family call.
    """
    if not os.path.isdir(sample_dir):
        return Counter(), "no dbCAN directory"
    cand = []
    for root, _dirs, files in os.walk(sample_dir):
        for f in files:
            if re.search(r"overview.*\.(tsv|txt|csv)$", f, re.I):
                cand.append(os.path.join(root, f))
    if not cand:
        return Counter(), "no overview table found"
    path = sorted(cand, key=len)[0]

    fams = Counter()
    n_rows = n_kept = 0
    with open(path) as fh:
        sniff = fh.readline()
        delim = "\t" if "\t" in sniff else ","
        fh.seek(0)
        rdr = csv.DictReader(fh, delimiter=delim)
        cols = rdr.fieldnames or []
        tool_col = next((c for c in cols if "tool" in c.lower()), None)
        fam_cols = [c for c in cols
                    if any(k in c.lower() for k in
                           ("hmm", "diamond", "dbcan_sub", "dbcansub", "ecami", "recommend"))]
        # Print exactly what was detected. dbCAN's overview column names have
        # changed between versions, and a silently mis-parsed table would put
        # wrong CAZyme counts into a manuscript - so make the parse auditable
        # rather than trusting it.
        print(f"     columns found : {cols}")
        print(f"     tool-count col: {tool_col or 'NONE (no >=2-tool filtering applied)'}")
        print(f"     family cols   : {fam_cols or 'NONE -> counts will be zero!'}")
        for row in rdr:
            n_rows += 1
            if tool_col:
                try:
                    if int(str(row[tool_col]).strip()) < 2:
                        continue
                except (ValueError, TypeError):
                    pass
            called = set()
            for c in fam_cols:
                v = (row.get(c) or "").strip()
                if not v or v in ("-", "N", "NA"):
                    continue
                for tok in re.split(r"[+|,;]", v):
                    m = re.match(r"([A-Z]{2,3})(\d+)", tok.strip())
                    if m:
                        called.add(m.group(1) + m.group(2))
            if called:
                n_kept += 1
            for f in called:
                fams[f] += 1
    print(f"     rows={n_rows}, proteins passing filter={n_kept}, distinct families={len(fams)}")
    return fams, os.path.relpath(path, sample_dir)


def find_faa(prokka_dir, sample):
    """Locate <sample>.faa under the Prokka output tree."""
    if not prokka_dir or not os.path.isdir(prokka_dir):
        return None
    for root, _d, files in os.walk(prokka_dir):
        for f in files:
            if f == f"{sample}.faa":
                return os.path.join(root, f)
    return None


def main():
    # .strip() because a mistyped shell line-continuation ("\" followed by a
    # space rather than a newline) silently prepends a space to the path,
    # which then fails to resolve in confusing ways.
    base = sys.argv[1].strip() if len(sys.argv) > 1 else "."
    # Optional 2nd arg: the Prokka output dir. IMPORTANT for a correct
    # annotation rate - the denominator must be the number of proteins Prokka
    # predicted, not a protein count from any other gene caller (e.g. the one
    # TYGS reports), or the percentage mixes two different gene-calling
    # pipelines and is not meaningful.
    prokka_dir = sys.argv[2].strip() if len(sys.argv) > 2 else None
    egg_base = os.path.join(base, "eggnog")
    caz_base = os.path.join(base, "dbcan")
    out = os.path.join(base, "summary")

    if not os.path.isdir(egg_base):
        sys.exit(f"No eggnog directory under {base} - is the path right?")
    # Only treat a subdirectory as a sample if it actually contains an
    # .emapper.annotations file - otherwise a stray folder (e.g. this
    # script's own "summary" output dir) gets picked up as a sample.
    samples = []
    for d in sorted(os.listdir(egg_base)):
        p = os.path.join(egg_base, d)
        if not os.path.isdir(p):
            continue
        if any(f.endswith(".emapper.annotations")
               for _r, _dd, fs in os.walk(p) for f in fs):
            samples.append(d)
    if not samples:
        sys.exit(
            f"No sample directories containing *.emapper.annotations were found under\n"
            f"  {egg_base}\n"
            f"Check the path. Pass the functional-annotation output directory as the\n"
            f"first argument (the one containing 'eggnog/' and 'dbcan/'), on ONE line:\n"
            f"  python3 17_summarise_functional.py <annot_dir> <prokka_dir>"
        )
    print(f"Samples: {', '.join(samples)}")
    os.makedirs(out, exist_ok=True)   # only once we know there is real input

    cog_counts = defaultdict(Counter)      # sample -> Counter(letter)
    kegg_counts = defaultdict(Counter)     # sample -> Counter(pathway)
    caz_from_egg = defaultdict(Counter)
    summary = {}

    for s in samples:
        ann = None
        for root, _d, files in os.walk(os.path.join(egg_base, s)):
            for f in files:
                if f.endswith(".emapper.annotations"):
                    ann = os.path.join(root, f)
        if not ann:
            print(f"  !! {s}: no .emapper.annotations found, skipping")
            continue

        n_rows = n_cog = n_ko = n_caz = 0
        for r in parse_emapper(ann):
            n_rows += 1
            cog = (r.get("COG_category") or "-").strip() or "-"
            if cog and cog != "-":
                n_cog += 1
                # a protein may carry several category letters; count each,
                # so column sums exceed the protein count (stated in the legend)
                for letter in cog:
                    if letter in COG_DESC:
                        cog_counts[s][letter] += 1
            else:
                cog_counts[s]["-"] += 1
            ko = (r.get("KEGG_ko") or "").strip()
            if ko and ko != "-":
                n_ko += 1
            paths = (r.get("KEGG_Pathway") or "").strip()
            if paths and paths != "-":
                for p in paths.split(","):
                    p = p.strip()
                    if p.startswith("ko"):
                        kegg_counts[s][p] += 1
            caz = (r.get("CAZy") or "").strip()
            if caz and caz != "-":
                n_caz += 1
                for tok in re.split(r"[,|]", caz):
                    m = re.match(r"([A-Z]{2,3})(\d+)", tok.strip())
                    if m:
                        caz_from_egg[s][m.group(1) + m.group(2)] += 1

        total_prot = count_faa(find_faa(prokka_dir, s))
        summary[s] = {
            "proteins_predicted": total_prot if total_prot else "",
            "annotation_rate_pct": (round(100.0 * n_rows / total_prot, 1)
                                    if total_prot else ""),
            "proteins_annotated": n_rows,
            "with_COG_category": n_cog,
            "with_KEGG_KO": n_ko,
            "with_CAZy_(eggNOG)": n_caz,
        }
        print(f"  {s}: {n_rows} annotated rows, {n_cog} with COG, {n_ko} with KEGG KO")

    # ---- dbCAN ----
    caz_counts = {}
    for s in samples:
        fams, note = parse_dbcan(os.path.join(caz_base, s))
        caz_counts[s] = fams
        summary.setdefault(s, {})["dbCAN_CAZymes"] = sum(fams.values())
        summary[s]["dbCAN_source"] = note
        print(f"  {s}: dbCAN -> {sum(fams.values())} CAZyme calls ({note})")
    if not any(caz_counts.values()):
        print("\n  !! dbCAN produced NO CAZyme calls for any sample.")
        print("     Check that script 13's dbCAN step actually ran and wrote output.")

    # ---------------------------------------------------------------- writers
    def write(name, header, rows):
        p = os.path.join(out, name)
        with open(p, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(rows)
        print(f"  wrote {p}")

    write("annotation_summary.csv",
          ["Sample", "Proteins predicted (Prokka)", "Proteins annotated (eggNOG)",
           "Annotation rate (%)", "With COG category", "With KEGG KO",
           "With CAZy (eggNOG)", "dbCAN CAZymes", "dbCAN source"],
          [[s,
            summary.get(s, {}).get("proteins_predicted", ""),
            summary.get(s, {}).get("proteins_annotated", 0),
            summary.get(s, {}).get("annotation_rate_pct", ""),
            summary.get(s, {}).get("with_COG_category", 0),
            summary.get(s, {}).get("with_KEGG_KO", 0),
            summary.get(s, {}).get("with_CAZy_(eggNOG)", 0),
            summary.get(s, {}).get("dbCAN_CAZymes", 0),
            summary.get(s, {}).get("dbCAN_source", "")] for s in samples])

    letters = [l for l in COG_DESC if any(cog_counts[s][l] for s in samples)]
    write("cog_counts.csv",
          ["COG category", "Description"] + samples,
          [[l, COG_DESC[l]] + [cog_counts[s][l] for s in samples] for l in letters])

    pct_rows = []
    for l in letters:
        row = [l, COG_DESC[l]]
        for s in samples:
            tot = sum(v for k, v in cog_counts[s].items() if k != "-")
            row.append(round(100.0 * cog_counts[s][l] / tot, 3) if tot else 0.0)
        pct_rows.append(row)
    write("cog_percent.csv", ["COG category", "Description"] + samples, pct_rows)

    allp = Counter()
    for s in samples:
        allp.update(kegg_counts[s])
    top = [p for p, _ in allp.most_common(60)]
    write("kegg_pathway_counts.csv", ["KEGG pathway"] + samples,
          [[p] + [kegg_counts[s][p] for s in samples] for p in top])

    src = caz_counts if any(caz_counts.values()) else caz_from_egg
    label = "dbCAN" if any(caz_counts.values()) else "eggNOG CAZy field (dbCAN empty)"
    print(f"  CAZyme tables built from: {label}")
    write("cazyme_class_counts.csv", ["CAZyme class"] + samples,
          [[c] + [sum(v for f, v in src[s].items() if re.match(rf"^{c}\d", f))
                  for s in samples] for c in CAZY_CLASSES])
    allf = Counter()
    for s in samples:
        allf.update(src[s])
    write("cazyme_family_counts.csv", ["CAZyme family"] + samples,
          [[f] + [src[s][f] for s in samples] for f, _ in allf.most_common()])

    print(f"\nDone. Send the contents of {out} back for figure generation.")


if __name__ == "__main__":
    main()
