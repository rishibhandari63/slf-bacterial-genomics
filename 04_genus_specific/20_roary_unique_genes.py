#!/usr/bin/env python3
"""
Find the genes that are present in one isolate and absent from every reference
genome in its genus-level Roary pan-genome, and profile them functionally.

This is the analysis that answers "what makes this novel species different from
its described congeners", as opposed to "how different is it".

Standard library only - runs on the Ceres login node.

Usage - ONE line:
    python3 20_roary_unique_genes.py <roary_dir> [emapper.annotations] [isolate_column]

e.g.
    python3 20_roary_unique_genes.py /project/kanglab/rishi.bhandari/agromyces/agromyces_refseq/core_phylo/roary_aln /project/kanglab/rishi.bhandari/funct_annotation_out/eggnog/Sample8/Sample8.emapper.annotations

If the isolate column cannot be auto-detected the script lists every genome
column it found, so you can pass the right name as the third argument.

IMPORTANT - locus tag matching. Roary's cells contain the locus tags assigned by
whichever Prokka run produced the GFFs you fed it. If your functional-annotation
Prokka run was SEPARATE from your Roary Prokka run, those tags will not match the
eggNOG query IDs and the functional join will find nothing. The script reports
the match rate explicitly and tells you what to do about it rather than silently
emitting an empty table.

Writes into <roary_dir>/unique_genes/:
    <isolate>_unique_genes.csv    - one row per isolate-specific cluster
    <isolate>_unique_locus_tags.txt - bare list, for pulling sequences
    <isolate>_absent_genes.csv    - clusters present in >=90% of references but
                                    absent from the isolate (candidate losses)
    <isolate>_cog_profile.csv     - COG profile of the unique genes (if joined)
"""
import csv
import os
import re
import sys
from collections import Counter

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

# Roary's fixed metadata columns; everything after these is a genome.
META = ["Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
        "Avg sequences per isolate", "Genome Fragment", "Order within Fragment",
        "Accessory Fragment", "Accessory Order with Fragment", "QC",
        "Min group size nuc", "Max group size nuc", "Avg group size nuc"]

COG_DESC = {
    "J": "Translation, ribosomal structure and biogenesis", "A": "RNA processing and modification",
    "K": "Transcription", "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics", "D": "Cell cycle control, cell division",
    "V": "Defense mechanisms", "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis", "N": "Cell motility", "Z": "Cytoskeleton",
    "W": "Extracellular structures", "U": "Intracellular trafficking, secretion",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "X": "Mobilome: prophages, transposons", "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism", "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism", "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism", "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only", "S": "Function unknown", "-": "Unassigned",
}


def load_eggnog(path):
    """query id -> {COG_category, Description, Preferred_name, KEGG_ko, KEGG_Pathway}"""
    if not path or not os.path.isfile(path):
        return {}
    out, header = {}, None
    with open(path, errors="replace") as fh:
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
            p = line.split("\t")
            if len(p) < len(header):
                p += [""] * (len(header) - len(p))
            d = dict(zip(header, p))
            out[d["query"]] = d
    return out


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    roary = sys.argv[1].strip()
    egg_path = sys.argv[2].strip() if len(sys.argv) > 2 else None
    forced_col = sys.argv[3].strip() if len(sys.argv) > 3 else None

    gpa = os.path.join(roary, "gene_presence_absence.csv")
    if not os.path.isfile(gpa):
        sys.exit(f"gene_presence_absence.csv not found in {roary}")

    with open(gpa, errors="replace") as fh:
        rdr = csv.DictReader(fh)
        fields = rdr.fieldnames or []
        rows = list(rdr)
    genomes = [f for f in fields if f not in META]
    print(f"Roary table: {len(rows)} gene clusters across {len(genomes)} genomes\n")

    # ---- identify the isolate column ----
    if forced_col:
        if forced_col not in genomes:
            sys.exit(f"Column {forced_col!r} not present. Genome columns are:\n  "
                     + "\n  ".join(genomes))
        iso = forced_col
    else:
        cands = [g for g in genomes if re.search(r"sample\s*\d+", g, re.I)]
        if len(cands) == 1:
            iso = cands[0]
        else:
            print("Could not uniquely auto-detect the isolate column.")
            print(f"Candidates matching 'Sample<N>': {cands or 'none'}")
            print("\nAll genome columns:")
            for g in genomes:
                print(f"  {g}")
            sys.exit("\nRerun passing the correct column name as the third argument.")
    refs = [g for g in genomes if g != iso]
    print(f"Isolate column : {iso}")
    print(f"Reference genomes: {len(refs)}\n")
    if len(refs) < 10:
        print(f"  !! Only {len(refs)} reference genomes. 'Unique' is a weak claim at this n -")
        print("     report the number of comparator genomes wherever you state it.\n")

    # ---- unique and absent clusters ----
    unique, absent = [], []
    for r in rows:
        in_iso = bool((r.get(iso) or "").strip())
        n_refs = sum(1 for g in refs if (r.get(g) or "").strip())
        if in_iso and n_refs == 0:
            unique.append(r)
        elif not in_iso and refs and n_refs / len(refs) >= 0.90:
            absent.append(r)

    print(f"Clusters unique to {iso}                    : {len(unique)}")
    print(f"Clusters in >=90% of references, absent here : {len(absent)}")

    # ---- functional join ----
    egg = load_eggnog(egg_path)
    matched = 0
    if egg:
        tags = []
        for r in unique:
            tags += [t.strip() for t in re.split(r"[\t,;]", r.get(iso, "")) if t.strip()]
        matched = sum(1 for t in tags if t in egg)
        rate = 100.0 * matched / len(tags) if tags else 0.0
        print(f"\neggNOG annotations loaded: {len(egg)} proteins")
        print(f"Locus-tag match rate     : {matched}/{len(tags)} ({rate:.1f}%)")
        if tags and rate < 50:
            print("\n  !! LOW MATCH RATE. The locus tags in the Roary table do not correspond")
            print("     to the eggNOG query IDs. This almost always means the Prokka run that")
            print("     produced your Roary GFFs was separate from the one that produced the")
            print("     proteins you annotated with eggNOG, so the two use different tags.")
            print("     Fix: take the .faa from the SAME Prokka run you gave Roary, pull the")
            print("     unique-gene proteins listed in the locus-tag file this script writes,")
            print("     and run eggNOG-mapper on just those (a few hundred proteins, minutes).")
            print("     The unique-gene list itself and its Prokka annotations are unaffected.")
    elif egg_path:
        print(f"\n  !! eggNOG file not readable: {egg_path}")

    out = os.path.join(roary, "unique_genes")
    os.makedirs(out, exist_ok=True)
    stem = re.sub(r"[^A-Za-z0-9]+", "_", iso).strip("_")

    def write(name, header, data):
        p = os.path.join(out, name)
        with open(p, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(data)
        print(f"  wrote {p}  ({len(data)} rows)")

    print()
    urows, all_tags = [], []
    cog = Counter()
    for r in unique:
        tags = [t.strip() for t in re.split(r"[\t,;]", r.get(iso, "")) if t.strip()]
        all_tags += tags
        e = next((egg[t] for t in tags if t in egg), {})
        c = (e.get("COG_category") or "").strip()
        for letter in (c if c and c != "-" else "-"):
            if letter in COG_DESC:
                cog[letter] += 1
        urows.append([
            r.get("Gene", ""), "; ".join(tags), r.get("Annotation", ""),
            e.get("Preferred_name", ""), c,
            "; ".join(COG_DESC.get(x, "") for x in c if x in COG_DESC),
            e.get("Description", ""), e.get("KEGG_ko", ""), e.get("KEGG_Pathway", ""),
        ])
    write(f"{stem}_unique_genes.csv",
          ["Roary cluster", "Locus tag(s)", "Prokka annotation", "Gene name",
           "COG category", "COG description", "eggNOG description",
           "KEGG KO", "KEGG pathway"], urows)

    p = os.path.join(out, f"{stem}_unique_locus_tags.txt")
    with open(p, "w") as fh:
        fh.write("\n".join(all_tags) + ("\n" if all_tags else ""))
    print(f"  wrote {p}  ({len(all_tags)} tags)")

    write(f"{stem}_absent_genes.csv",
          ["Roary cluster", "Prokka annotation", "No. isolates with gene"],
          [[r.get("Gene", ""), r.get("Annotation", ""), r.get("No. isolates", "")]
           for r in absent])

    if cog:
        tot = sum(v for k, v in cog.items() if k != "-")
        write(f"{stem}_cog_profile.csv",
              ["COG category", "Description", "Unique genes", "% of COG-assigned"],
              [[k, COG_DESC[k], v, round(100.0 * v / tot, 1) if tot and k != "-" else ""]
               for k, v in cog.most_common()])

    # ---- how much is hypothetical? the honest headline ----
    hyp = sum(1 for r in unique
              if re.search(r"hypothetical|unknown|uncharacter", r.get("Annotation", ""), re.I))
    if unique:
        print(f"\n  Of {len(unique)} isolate-specific clusters, {hyp} "
              f"({100.0 * hyp / len(unique):.0f}%) are annotated only as hypothetical/unknown.")
        print("  Report this fraction - it is the honest measure of how much of the")
        print("  novelty is interpretable versus merely detected.")

    print(f"\nDone. Send the contents of {out} back for the write-up.")
    print("\nREMINDER: absence in a draft reference assembly can be an assembly gap.")
    print("Any gene you highlight in the manuscript should be confirmed absent by")
    print("tblastn against the reference ASSEMBLIES, not just missing from annotations.")


if __name__ == "__main__":
    main()
