#!/usr/bin/env python3
"""Compute POCP and AAI from reciprocal DIAMOND blastp output.

POCP (percentage of conserved proteins), Qin et al. 2014, J Bacteriol 196:2210.
    A protein counts as conserved when its best hit has
        E-value  < 1e-5
        identity > 40%
        alignment length > 50% of the QUERY protein length
    POCP = (C1 + C2) / (T1 + T2) * 100
    where C1 = conserved proteins of genome 1 against genome 2,
          T1 = total proteins of genome 1.
    Two genomes are conventionally regarded as belonging to the same genus
    when POCP > 50%.

AAI (average amino-acid identity) is computed here from one-way best hits
    passing identity > 30% and query coverage > 50%, then averaged over both
    directions. Reported for corroboration only; the genus boundary is usually
    quoted around 65%, but it is softer than POCP and varies between lineages.

Standard library only.
"""
import argparse
import csv
import glob
import os
import sys


def count_proteins(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def parse_hits(path):
    """Return {query: (pident, length, qlen)} keeping the best bitscore per query."""
    best = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            q = f[0]
            try:
                pident, length, qlen = float(f[2]), int(f[3]), int(f[4])
                bits = float(f[7])
            except ValueError:
                continue
            prev = best.get(q)
            if prev is None or bits > prev[3]:
                best[q] = (pident, length, qlen, bits)
    return best


def conserved_count(hits):
    """POCP criteria: identity > 40% and alignment length > 50% of query length."""
    n = 0
    for pident, length, qlen, _ in hits.values():
        if qlen and pident > 40.0 and length > 0.5 * qlen:
            n += 1
    return n


def aai_values(hits):
    """AAI criteria: identity > 30% and query coverage > 50%."""
    vals = []
    for pident, length, qlen, _ in hits.values():
        if qlen and pident > 30.0 and length > 0.5 * qlen:
            vals.append(pident)
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--faa-dir", required=True)
    ap.add_argument("--hits-dir", required=True)
    ap.add_argument("--focal", required=True,
                    help="basename of the isolate proteome, e.g. Sample4")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    faas = sorted(glob.glob(os.path.join(args.faa_dir, "*.faa")))
    if not faas:
        sys.exit(f"ERROR: no .faa files in {args.faa_dir}")

    totals = {}
    for f in faas:
        b = os.path.basename(f)[:-4]
        totals[b] = count_proteins(f)
    if args.focal not in totals:
        sys.exit(f"ERROR: focal proteome '{args.focal}' not among {sorted(totals)}")

    print(f"Proteomes ({len(totals)}):")
    for b in sorted(totals):
        mark = "  <- focal" if b == args.focal else ""
        print(f"  {b:32s} {totals[b]:6d} proteins{mark}")

    rows = []
    others = [b for b in sorted(totals) if b != args.focal]
    missing = []

    for other in others:
        fwd_p = os.path.join(args.hits_dir, f"{args.focal}__vs__{other}.tsv")
        rev_p = os.path.join(args.hits_dir, f"{other}__vs__{args.focal}.tsv")
        if not (os.path.exists(fwd_p) and os.path.exists(rev_p)):
            missing.append(other)
            continue

        fwd, rev = parse_hits(fwd_p), parse_hits(rev_p)
        c1, c2 = conserved_count(fwd), conserved_count(rev)
        t1, t2 = totals[args.focal], totals[other]
        pocp = (c1 + c2) / (t1 + t2) * 100.0 if (t1 + t2) else 0.0

        v = aai_values(fwd) + aai_values(rev)
        aai = sum(v) / len(v) if v else 0.0

        rows.append({
            "Comparator": other,
            "Comparator proteins": t2,
            "Conserved (focal->comp)": c1,
            "Conserved (comp->focal)": c2,
            "POCP (%)": round(pocp, 2),
            "Same genus by POCP": "yes" if pocp > 50 else "no",
            "AAI (%)": round(aai, 2),
            "AAI pairs": len(v),
        })

    if missing:
        print("\nWARNING: no DIAMOND output for: " + ", ".join(missing), file=sys.stderr)
        print("These comparators are absent from the table. Check step 3 of the "
              "slurm script completed for them.", file=sys.stderr)

    if not rows:
        sys.exit("ERROR: no comparisons could be computed - nothing written.")

    rows.sort(key=lambda r: -r["POCP (%)"])
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"\nPOCP / AAI of {args.focal} against each type strain")
    print(f"  {'comparator':32s} {'POCP':>8s} {'genus?':>8s} {'AAI':>8s}")
    for r in rows:
        print(f"  {r['Comparator']:32s} {r['POCP (%)']:8.2f} "
              f"{r['Same genus by POCP']:>8s} {r['AAI (%)']:8.2f}")

    print(f"\nWritten: {args.out}")
    print("\nInterpretation: POCP > 50% indicates the same genus. If Sample4 "
          "exceeds 50% against Massilia timonae AND Telluria mixta, that is "
          "the expected result - it confirms the two genera are not separable "
          "by this criterion, which is exactly why the nomenclature is "
          "contested. What matters then is that Sample4 does NOT fall closer "
          "to Mokoshia (M. eurypsychrophila) or Zemynaea (M. arenosa).")


if __name__ == "__main__":
    main()
