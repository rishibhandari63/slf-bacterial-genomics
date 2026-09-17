#!/usr/bin/env python3
"""
Summarise the MGE / defense-system outputs (script 15) into tidy CSV tables.

Standard library only - runs on the Ceres login node, no conda env needed.

Parses, per sample:
  geNomad       - predicted viruses (prophages) and plasmid contigs
  MOB-suite     - contig_report.txt: chromosome vs plasmid assignment, clusters, mobility
  AMRFinderPlus - amrfinder.tsv: hits split by element type (AMR / STRESS / VIRULENCE)
  DefenseFinder - defense_finder_systems.tsv: anti-phage systems by type/subtype

Each parser prints the file it used and the columns it found, so the parse can
be verified rather than trusted (tool output formats drift between versions).

Usage - ONE line:
    python3 19_summarise_mge_defense.py <mge_defense_out_dir>

e.g.
    python3 19_summarise_mge_defense.py /project/kanglab/rishi.bhandari/mge_defense_out

Writes into <dir>/summary/:
    mge_summary.csv          - one row per genome, all headline counts
    defense_system_counts.csv- defense system type x genome matrix
    defense_systems.csv      - every detected system, one row each
    amr_hits.csv             - every AMRFinderPlus hit, one row each
    mobile_elements.csv      - every predicted prophage / plasmid, one row each
"""
import csv
import os
import re
import sys
from collections import Counter, defaultdict


def find_file(root, *patterns):
    """First file under root matching any regex, shallowest path first."""
    hits = []
    if not os.path.isdir(root):
        return None
    for dirpath, _d, files in os.walk(root):
        for f in files:
            for p in patterns:
                if re.search(p, f, re.I):
                    hits.append(os.path.join(dirpath, f))
    return sorted(hits, key=lambda p: (p.count(os.sep), len(p)))[0] if hits else None


def read_table(path):
    """Read a TSV/CSV into list-of-dicts; returns ([], []) if unreadable."""
    if not path or not os.path.isfile(path):
        return [], []
    with open(path, errors="replace") as fh:
        first = fh.readline()
        if not first.strip():
            return [], []
        delim = "\t" if "\t" in first else ","
        fh.seek(0)
        rdr = csv.DictReader(fh, delimiter=delim)
        rows = [r for r in rdr]
        return rows, (rdr.fieldnames or [])


def col(fields, *cands):
    """Find a column by case-insensitive substring match."""
    for c in cands:
        for f in fields:
            if c.lower() in f.lower():
                return f
    return None


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python3 19_summarise_mge_defense.py <mge_defense_out_dir>")
    base = sys.argv[1].strip()
    if not os.path.isdir(base):
        sys.exit(f"Not a directory: {base}")

    samples = sorted(d for d in os.listdir(base)
                     if os.path.isdir(os.path.join(base, d)) and d != "summary")
    if not samples:
        sys.exit(f"No sample directories under {base}")
    print(f"Samples: {', '.join(samples)}\n")

    summary = {}
    defense_rows, amr_rows, mge_rows = [], [], []
    defense_counts = defaultdict(Counter)
    printed_cols = set()

    for s in samples:
        sd = os.path.join(base, s)
        rec = {}

        # ---------------- geNomad ----------------
        vpath = find_file(os.path.join(sd, "genomad"), r"virus_summary\.tsv$")
        ppath = find_file(os.path.join(sd, "genomad"), r"plasmid_summary\.tsv$")
        vrows, vf = read_table(vpath)
        prows, pf = read_table(ppath)
        rec["genomad_viruses"] = len(vrows)
        rec["genomad_plasmids"] = len(prows)
        if vpath and "genomad" not in printed_cols:
            print(f"  [geNomad] virus table: {os.path.basename(vpath)}  cols={vf[:6]}")
            printed_cols.add("genomad")
        for r in vrows:
            mge_rows.append({
                "Sample": s, "Element": "prophage/virus",
                "Sequence": r.get(col(vf, "seq_name") or "", ""),
                "Length (bp)": r.get(col(vf, "length") or "", ""),
                "Score": r.get(col(vf, "virus_score", "score") or "", ""),
                "Detail": r.get(col(vf, "taxonomy") or "", ""),
            })
        for r in prows:
            mge_rows.append({
                "Sample": s, "Element": "plasmid (geNomad)",
                "Sequence": r.get(col(pf, "seq_name") or "", ""),
                "Length (bp)": r.get(col(pf, "length") or "", ""),
                "Score": r.get(col(pf, "plasmid_score", "score") or "", ""),
                "Detail": f"conj_genes={r.get(col(pf, 'conjugation') or '', '')}",
            })

        # ---------------- MOB-suite ----------------
        cpath = find_file(os.path.join(sd, "mob_recon"), r"contig_report\.(txt|tsv)$")
        crows, cf = read_table(cpath)
        mol = col(cf, "molecule_type")
        clus = col(cf, "primary_cluster_id")
        mob = col(cf, "predicted_mobility")
        plas = [r for r in crows if mol and str(r.get(mol, "")).lower().startswith("plasmid")]
        rec["mob_plasmid_contigs"] = len(plas)
        rec["mob_plasmid_clusters"] = len({r.get(clus) for r in plas if r.get(clus)}) if clus else ""
        rec["mob_chromosome_contigs"] = sum(
            1 for r in crows if mol and str(r.get(mol, "")).lower().startswith("chromosome"))
        if cpath and "mob" not in printed_cols:
            print(f"  [MOB-suite] {os.path.basename(cpath)}  molecule_type col={mol!r}, "
                  f"cluster col={clus!r}")
            printed_cols.add("mob")
        for r in plas:
            mge_rows.append({
                "Sample": s, "Element": "plasmid (MOB-suite)",
                "Sequence": r.get(col(cf, "contig_id") or "", ""),
                "Length (bp)": r.get(col(cf, "size") or "", ""),
                "Score": r.get(clus or "", ""),
                "Detail": f"mobility={r.get(mob or '', '')}",
            })

        # ---------------- AMRFinderPlus ----------------
        apath = os.path.join(sd, "amrfinder.tsv")
        arows, af = read_table(apath)
        et = col(af, "Element type")
        types = Counter(r.get(et, "?") for r in arows) if et else Counter()
        rec["amr_total"] = len(arows)
        rec["amr_AMR"] = types.get("AMR", 0)
        rec["amr_STRESS"] = types.get("STRESS", 0)
        rec["amr_VIRULENCE"] = types.get("VIRULENCE", 0)
        if arows and "amr" not in printed_cols:
            print(f"  [AMRFinder] element-type col={et!r}")
            printed_cols.add("amr")
        for r in arows:
            amr_rows.append({
                "Sample": s,
                # Contig id matters: it is what tells you whether a resistance
                # gene sits on a plasmid contig or on the chromosome, which is
                # the difference between "carries resistance" and "carries
                # MOBILE resistance". Cross-reference against the plasmid
                # contigs listed in mobile_elements.csv.
                "Contig id": r.get(col(af, "Contig id", "Contig") or "", ""),
                "Start": r.get(col(af, "^Start$", "Start") or "", ""),
                "Gene symbol": r.get(col(af, "Gene symbol") or "", ""),
                "Sequence name": r.get(col(af, "Sequence name") or "", ""),
                "Element type": r.get(et or "", ""),
                "Element subtype": r.get(col(af, "Element subtype") or "", ""),
                "Class": r.get(col(af, "^Class") or col(af, "Class") or "", ""),
                "Subclass": r.get(col(af, "Subclass") or "", ""),
                "% identity": r.get(col(af, "identity") or "", ""),
            })

        # ---------------- DefenseFinder ----------------
        dpath = find_file(os.path.join(sd, "defensefinder"), r"systems\.tsv$")
        drows, df = read_table(dpath)
        tcol = col(df, "^type$", "type")
        scol = col(df, "subtype")
        rec["defense_systems"] = len(drows)
        rec["defense_types"] = len({r.get(tcol) for r in drows if tcol and r.get(tcol)})
        if dpath and "df" not in printed_cols:
            print(f"  [DefenseFinder] {os.path.basename(dpath)}  type col={tcol!r}, "
                  f"subtype col={scol!r}")
            printed_cols.add("df")
        for r in drows:
            t = r.get(tcol, "?") if tcol else "?"
            defense_counts[s][t] += 1
            defense_rows.append({
                "Sample": s,
                "System type": t,
                "Subtype": r.get(scol, "") if scol else "",
                "Genes": r.get(col(df, "genes_count", "gene_count") or "", ""),
                "System id": r.get(col(df, "sys_id") or "", ""),
            })

        summary[s] = rec
        print(f"  {s}: {rec['genomad_viruses']} prophage/virus, "
              f"{rec['genomad_plasmids']} geNomad plasmid, "
              f"{rec['mob_plasmid_contigs']} MOB plasmid contigs, "
              f"{rec['amr_total']} AMR/stress/virulence hits, "
              f"{rec['defense_systems']} defense systems")

    out = os.path.join(base, "summary")
    os.makedirs(out, exist_ok=True)

    def write(name, header, rows):
        p = os.path.join(out, name)
        with open(p, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(rows)
        print(f"  wrote {p}  ({len(rows)} rows)")

    print()
    keys = ["genomad_viruses", "genomad_plasmids", "mob_plasmid_contigs",
            "mob_plasmid_clusters", "mob_chromosome_contigs",
            "amr_total", "amr_AMR", "amr_STRESS", "amr_VIRULENCE",
            "defense_systems", "defense_types"]
    write("mge_summary.csv", ["Sample"] + keys,
          [[s] + [summary[s].get(k, "") for k in keys] for s in samples])

    all_types = sorted({t for s in samples for t in defense_counts[s]})
    write("defense_system_counts.csv", ["Defense system type"] + samples,
          [[t] + [defense_counts[s][t] for s in samples] for t in all_types])

    for name, rows in (("defense_systems.csv", defense_rows),
                       ("amr_hits.csv", amr_rows),
                       ("mobile_elements.csv", mge_rows)):
        if rows:
            write(name, list(rows[0].keys()), [list(r.values()) for r in rows])
        else:
            write(name, ["(no rows detected)"], [])

    print(f"\nDone. Send the contents of {out} back for figures and write-up.")


if __name__ == "__main__":
    main()
