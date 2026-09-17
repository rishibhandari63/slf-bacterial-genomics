#!/usr/bin/env python3
"""
Summarise antiSMASH output into tidy CSV tables.

Standard library only - runs on the Ceres login node with no conda env.

antiSMASH's output layout and field names shift between major versions, so
this script reads the per-region GenBank files (the most stable artefact),
tries several sources for the MIBiG/knownclusterblast similarity, and PRINTS
which source it used plus a sample of what it parsed. Verify that diagnostic
output before trusting the numbers.

Usage - type it on ONE line:
    python3 18_summarise_antismash.py <antismash_out_dir>

e.g.
    python3 18_summarise_antismash.py /project/kanglab/rishi.bhandari/antismash_out

Writes into <antismash_out_dir>/summary/:
    bgc_regions.csv   - one row per BGC region (sample, region, product, length, contig edge, best MIBiG hit)
    bgc_type_counts.csv - BGC product type x sample count matrix
    bgc_summary.csv   - per sample: total regions, complete vs contig-edge, distinct types
"""
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict


def parse_region_gbk(path):
    """Extract region-level info from one .regionNNN.gbk file."""
    txt = open(path, errors="replace").read()

    # total length from the LOCUS line
    m = re.search(r"^LOCUS\s+\S+\s+(\d+)\s+bp", txt, re.M)
    length = int(m.group(1)) if m else None

    # the 'region' feature block: from '     region ' to the next feature at
    # the same indent level
    rm = re.search(r"^ {5}region\s+.*?\n(.*?)(?=^ {5}\S)", txt, re.M | re.S)
    block = rm.group(1) if rm else ""

    def quals(name):
        return re.findall(rf'/{name}="([^"]*)"', block, re.S)

    products = [p.strip() for p in quals("product")]
    edge_q = quals("contig_edge")
    edge = edge_q[0] if edge_q else "__ABSENT__"
    num = (quals("region_number") or ["?"])[0]

    # knownclusterblast qualifiers, when present, look like:
    #   /knownclusterblast="1. BGC0000001_c1  Abyssomicin (55% of genes ...)"
    kcb = quals("knownclusterblast")
    n_cds = len(re.findall(r"^ {5}CDS\s+", txt, re.M))
    return {
        "n_cds": n_cds,
        "length": length,
        "products": products,
        "contig_edge": edge,
        "region_number": num,
        "kcb_qualifiers": kcb,
        "raw_block_sample": block[:400],
    }


def best_hit_from_kcb_quals(kcb):
    """Pick the top-ranked knownclusterblast qualifier and its similarity."""
    if not kcb:
        return "", ""
    first = sorted(kcb, key=lambda s: int(re.match(r"\s*(\d+)", s).group(1))
                   if re.match(r"\s*(\d+)", s) else 999)[0]
    pct = ""
    pm = re.search(r"(\d+)%", first)
    if pm:
        pct = pm.group(1)
    name = re.sub(r"^\s*\d+\.\s*", "", first)
    name = re.sub(r"\s*\(\d+% of genes.*?\)\s*$", "", name).strip()
    # the qualifier joins accession and description with a tab
    parts = [p.strip() for p in re.split(r"\t+", name) if p.strip()]
    if len(parts) >= 2:
        return f"{parts[1]} ({parts[0]})", pct
    return name, pct


def best_hit_from_json(sample_dir, record_id, region_no):
    """
    Fall back to the antiSMASH JSON, which carries knownclusterblast results
    under each record's modules. Structure varies by version, so this walks
    defensively and returns ('', '') rather than raising.
    """
    cands = [os.path.join(sample_dir, f) for f in os.listdir(sample_dir)
             if f.endswith(".json")]
    if not cands:
        return "", "", "no json"
    try:
        data = json.load(open(sorted(cands, key=len)[0], errors="replace"))
    except Exception as e:
        return "", "", f"json unreadable ({e.__class__.__name__})"
    try:
        for rec in data.get("records", []):
            mods = rec.get("modules", {})
            kcb = mods.get("antismash.modules.clusterblast", {})
            known = kcb.get("knowncluster", {})
            results = known.get("results", [])
            for res in results:
                if str(res.get("region_number")) != str(region_no):
                    continue
                ranking = res.get("ranking") or []
                if not ranking:
                    return "", "", "json ok (no hits)"
                top = ranking[0]
                desc = top[0] if isinstance(top, list) else {}
                score = top[1] if isinstance(top, list) and len(top) > 1 else {}
                name = (desc.get("description") or desc.get("accession") or "")
                hits = score.get("hits") if isinstance(score, dict) else None
                return name, (str(hits) if hits is not None else ""), "json ok"
    except Exception as e:
        return "", "", f"json walk failed ({e.__class__.__name__})"
    return "", "", "json ok (region not found)"


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python3 18_summarise_antismash.py <antismash_out_dir>")
    base = sys.argv[1].strip()
    if not os.path.isdir(base):
        sys.exit(f"Not a directory: {base}")

    samples = []
    for d in sorted(os.listdir(base)):
        p = os.path.join(base, d)
        if not os.path.isdir(p) or d == "summary":
            continue
        if any(re.search(r"\.region\d+\.gbk$", f) for f in os.listdir(p)):
            samples.append(d)
        else:
            print(f"  note: {d} has no .region*.gbk files "
                  f"(contains: {sorted(os.listdir(p))[:6]}...)")
    if not samples:
        sys.exit(
            f"No sample directories with .region*.gbk found under {base}.\n"
            "antiSMASH writes one directory per sample; check the path."
        )
    print(f"Samples: {', '.join(samples)}\n")

    rows = []
    type_counts = defaultdict(Counter)
    per_sample = {}
    kcb_source_note = set()

    for s in samples:
        sdir = os.path.join(base, s)
        regs = sorted(f for f in os.listdir(sdir) if re.search(r"\.region\d+\.gbk$", f))
        n_edge = 0
        n_absent = 0
        for f in regs:
            info = parse_region_gbk(os.path.join(sdir, f))
            prod = "; ".join(info["products"]) or "unknown"
            raw_edge = info["contig_edge"]
            if raw_edge == "__ABSENT__":
                edge, edge_label = False, "not recorded"
                n_absent += 1
            else:
                edge = raw_edge.lower() in ("true", "yes", "1")
                edge_label = "yes" if edge else "no"
            if edge:
                n_edge += 1

            name, pct = best_hit_from_kcb_quals(info["kcb_qualifiers"])
            metric = ""
            if name:
                metric = "percent_of_genes"
                kcb_source_note.add("region .gbk /knownclusterblast qualifier (a true %)")
            else:
                rec_id = f.split(".region")[0]
                name, pct, note = best_hit_from_json(sdir, rec_id, info["region_number"])
                if name:
                    # NOTE: this is the JSON "hits" field = NUMBER OF GENES with a
                    # BLAST hit to the MIBiG cluster. It is NOT a percentage. An
                    # earlier version of this script mislabelled it as one.
                    metric = "gene_hit_count"
                kcb_source_note.add(f"antiSMASH JSON hit COUNT (not a %) [{note}]")

            rows.append({
                "Sample": s,
                "Region file": f,
                "Region": info["region_number"],
                "Product type(s)": prod,
                "Length (bp)": info["length"] or "",
                "CDS in region": info["n_cds"],
                "On contig edge": edge_label,
                "Closest MIBiG cluster": name,
                "MIBiG match value": pct,
                "MIBiG match metric": metric,
                "Genes hit / CDS in region (%)": (
                    round(100.0 * int(pct) / info["n_cds"], 1)
                    if metric == "gene_hit_count" and pct and info["n_cds"] else ""),
            })
            for p in (info["products"] or ["unknown"]):
                type_counts[s][p] += 1

        per_sample[s] = {
            "regions": len(regs),
            "contig_edge": n_edge,
            "complete": len(regs) - n_edge,
            "distinct_types": len(type_counts[s]),
        }
        msg = (f"  {s}: {len(regs)} BGC regions "
               f"({len(regs) - n_edge} complete, {n_edge} on contig edge), "
               f"{len(type_counts[s])} distinct product types")
        if n_absent:
            msg += f"  [!! {n_absent} region(s) had NO contig_edge qualifier - 'complete' is unverified]"
        print(msg)

    print("\n  MIBiG similarity read from: " + "; ".join(sorted(kcb_source_note)))
    if rows:
        print("  example parsed region:")
        for k, v in rows[0].items():
            print(f"     {k}: {v}")

    out = os.path.join(base, "summary")
    os.makedirs(out, exist_ok=True)

    def write(name, header, data):
        p = os.path.join(out, name)
        with open(p, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(data)
        print(f"  wrote {p}")

    print()
    write("bgc_regions.csv", list(rows[0].keys()), [list(r.values()) for r in rows])

    all_types = sorted({t for s in samples for t in type_counts[s]})
    write("bgc_type_counts.csv", ["BGC product type"] + samples,
          [[t] + [type_counts[s][t] for s in samples] for t in all_types])

    write("bgc_summary.csv",
          ["Sample", "Total BGC regions", "Complete", "On contig edge", "Distinct product types"],
          [[s, per_sample[s]["regions"], per_sample[s]["complete"],
            per_sample[s]["contig_edge"], per_sample[s]["distinct_types"]] for s in samples])

    print(f"\nDone. Send the contents of {out} back for figures and write-up.")
    print("NOTE: BGCs flagged 'On contig edge' are truncated by the assembly, so")
    print("      counts from draft genomes are lower bounds - report that caveat.")


if __name__ == "__main__":
    main()
