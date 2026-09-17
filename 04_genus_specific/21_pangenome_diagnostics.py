#!/usr/bin/env python3
"""
Diagnose a collapsed Roary core genome.

A core of tens rather than ~1000 genes has two common causes, and they need
different fixes:

  (A) OVER-SPLITTING. Roary's default 95% BLASTP identity is too strict across
      species boundaries within a genus, so true orthologs land in separate
      clusters. Signature: the core stays small even when you relax the
      presence threshold, and most clusters contain very few genomes.
      Fix: rerun at lower identity (roary -i 80) or use Panaroo.

  (B) A FEW BAD GENOMES. The core requires presence in 99% of genomes, so with
      153 genomes a single fragmented or misassigned assembly missing a gene
      removes it from the core. Signature: the core recovers sharply as the
      threshold drops, and a handful of genomes carry far fewer clusters than
      the rest. Fix: QC-filter the reference set, then rerun.

Both can apply at once. This script quantifies each without recomputing
anything - it only reads gene_presence_absence.csv.

Standard library only. Usage - ONE line:
    python3 21_pangenome_diagnostics.py <roary_dir>

Writes <roary_dir>/diagnostics/:
    core_vs_threshold.csv   - core size at a range of presence thresholds
    genes_per_genome.csv    - clusters per genome, flagging low outliers
"""
import csv
import os
import sys
from collections import Counter

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

META = {"Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
        "Avg sequences per isolate", "Genome Fragment", "Order within Fragment",
        "Accessory Fragment", "Accessory Order with Fragment", "QC",
        "Min group size nuc", "Max group size nuc", "Avg group size nuc"}

THRESHOLDS = [0.99, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.30, 0.15]


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python3 21_pangenome_diagnostics.py <roary_dir>")
    roary = sys.argv[1].strip()
    gpa = os.path.join(roary, "gene_presence_absence.csv")
    if not os.path.isfile(gpa):
        sys.exit(f"gene_presence_absence.csv not found in {roary}")

    print("Reading (streaming - large tables take a minute)...")
    with open(gpa, errors="replace") as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        gidx = [i for i, h in enumerate(header) if h not in META]
        genomes = [header[i] for i in gidx]
        n_gen = len(genomes)

        per_genome = Counter()          # genome -> clusters it appears in
        presence_hist = Counter()       # n_genomes_in_cluster -> cluster count
        n_clusters = 0
        for row in rdr:
            if not row:
                continue
            n_clusters += 1
            present = 0
            for i in gidx:
                if i < len(row) and row[i].strip():
                    per_genome[header[i]] += 1
                    present += 1
            presence_hist[present] += 1

    print(f"\nGenomes : {n_gen}")
    print(f"Clusters: {n_clusters:,}")

    # ---- core size as a function of presence threshold ----
    print("\nCore size vs presence threshold")
    print("  threshold   genomes required   core clusters")
    rows_out = []
    for t in THRESHOLDS:
        need = int(round(t * n_gen))
        core = sum(c for k, c in presence_hist.items() if k >= need)
        rows_out.append([f"{t:.0%}", need, core])
        print(f"  {t:>7.0%}   {need:>16}   {core:>13,}")

    core99 = rows_out[0][2]
    core50 = next(r[2] for r in rows_out if r[0] == "50%")
    core15 = rows_out[-1][2]

    # ---- clusters per genome: find the stragglers ----
    counts = sorted(per_genome.items(), key=lambda kv: kv[1])
    vals = sorted(per_genome.values())
    med = vals[len(vals) // 2] if vals else 0
    low = [(g, c) for g, c in counts if med and c < 0.75 * med]

    print(f"\nClusters per genome: median {med:,}, "
          f"min {vals[0]:,}, max {vals[-1]:,}" if vals else "")
    if low:
        print(f"\n  {len(low)} genome(s) carry <75% of the median cluster count.")
        print("  These are the most likely candidates for fragmented, incomplete or")
        print("  misassigned assemblies dragging the core down:")
        for g, c in low[:15]:
            print(f"     {c:>7,}  {g}")
        if len(low) > 15:
            print(f"     ... and {len(low) - 15} more (see genes_per_genome.csv)")
    else:
        print("\n  No genome carries markedly fewer clusters than the others.")

    # ---- verdict ----
    print("\n" + "=" * 68)
    print("INTERPRETATION")
    print("=" * 68)
    recovers = core50 >= 8 * max(core99, 1)
    singleton_frac = 100.0 * sum(c for k, c in presence_hist.items() if k <= 2) / max(n_clusters, 1)
    print(f"  Core at 99% presence            : {core99:,}")
    print(f"  Core at 50% presence            : {core50:,}")
    print(f"  Clusters in <=2 genomes         : {singleton_frac:.1f}% of all clusters")
    print()
    if singleton_frac > 60:
        print("  (A) OVER-SPLITTING is present. A large majority of clusters contain")
        print("      one or two genomes, which is the signature of an identity")
        print("      threshold too strict for between-species comparison.")
        print("      -> Rerun the pan-genome with `roary -i 80` (or lower), or with")
        print("         Panaroo, which is designed for diverse gene sets.")
    if low:
        print("  (B) LOW-QUALITY OR MISASSIGNED GENOMES are also contributing: the")
        print("      genomes listed above carry far fewer clusters than the median.")
        print("      -> QC-filter the reference set (CheckM2 completeness >=95%,")
        print("         contamination <=5%) and confirm each is really in this genus")
        print("         (GTDB-Tk), then rerun.")
    if recovers and not low:
        print("  The core recovers sharply as the threshold relaxes, which points to")
        print("  a strict-threshold effect rather than clustering failure alone.")
        print("  -> Consider `roary -cd 95` in addition to lowering identity.")
    if not low and singleton_frac <= 60:
        print("  Neither signature is strong; inspect genes_per_genome.csv by hand.")
    print()
    print("  NOTE: whatever you change, the core alignment feeding your RAxML tree")
    print("  is produced by the same clustering. A core of tens of genes is too thin")
    print("  a basis for a genus-level phylogeny - rerun the tree from the corrected")
    print("  pan-genome as well, not just the uniqueness analysis.")

    out = os.path.join(roary, "diagnostics")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "core_vs_threshold.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Presence threshold", "Genomes required", "Core clusters"])
        w.writerows(rows_out)
    with open(os.path.join(out, "genes_per_genome.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Genome", "Clusters present", "Fraction of median"])
        for g, c in counts:
            w.writerow([g, c, round(c / med, 3) if med else ""])
    print(f"\n  wrote {out}/core_vs_threshold.csv")
    print(f"  wrote {out}/genes_per_genome.csv")


if __name__ == "__main__":
    main()
