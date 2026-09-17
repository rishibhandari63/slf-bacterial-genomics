#!/bin/bash
#SBATCH --job-name=unique_genes
#SBATCH -A kanglab
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -o unique_genes_%j.out
#SBATCH -e unique_genes_%j.err

set -euo pipefail

# Identify the proteins an isolate carries that have NO counterpart in any
# reference genome of its genus, by direct homology search rather than by
# pan-genome clustering.
#
# WHY NOT ROARY HERE: Roary's clustering collapsed on these datasets (Agromyces
# core = 26 genes across 153 genomes) because a fixed 95% identity threshold
# cannot group orthologs across species boundaries within a genus. It reported
# 51% of Sample8's proteome as "unique"; this approach reports 0.8%.
#
# TO SWITCH BETWEEN ISOLATES: change the SAMPLE line below. Nothing else.

# ============================================================
# Which isolate to run. Default is Sample8; override at submit time
# WITHOUT editing this file:
#     sbatch 22_unique_genes_diamond.slurm                      # Sample8
#     sbatch --export=ALL,SAMPLE=Sample4 22_unique_genes_diamond.slurm
# ============================================================
SAMPLE="${SAMPLE:-Sample8}"

PROJECT=/project/kanglab/rishi.bhandari

# The isolate proteome MUST be the one eggNOG-mapper annotated, or the
# functional join silently returns nothing: Prokka assigns fresh locus tags on
# every run, so the .faa from the Roary Prokka run has different identifiers
# from the .faa under PROKKA_DIR even though the gene calls are identical.
ISOLATE_FAA=$PROJECT/PROKKA_DIR/prokka_out/$SAMPLE/$SAMPLE.faa
EGGNOG=$PROJECT/funct_annotation_out/eggnog/$SAMPLE/$SAMPLE.emapper.annotations

# Genus reference proteomes - searched recursively for *.faa
case "$SAMPLE" in
  Sample8) REF_DIR=$PROJECT/agromyces/agromyces_refseq ;;
  # Sample4's comparators are MASSILIA, not Telluria. Point this at whichever
  # directory holds the proteomes for the same genome set used for the core
  # phylogeny, so the tree and the gene-content analysis describe one genus.
  Sample4) REF_DIR=$PROJECT/massilia/faa ;;
  *) echo "ERROR: no REF_DIR configured for $SAMPLE - add a case below." >&2; exit 1 ;;
esac

OUTDIR=$PROJECT/unique_genes_$SAMPLE
DIAMOND=$PROJECT/envs/funct_env/bin/diamond

# A protein is ABSENT from the genus if it has no hit at this E-value;
# DIVERGENT if its best hit falls below this identity. Both are reported.
EVALUE=1e-5
DIVERGENT_MAX_ID=40

# -----------------------------------------------------------------
mkdir -p "$OUTDIR"
cd "$OUTDIR"

for f in "$ISOLATE_FAA" "$DIAMOND"; do
  [ -s "$f" ] || { echo "ERROR: not found: $f" >&2; exit 1; }
done
[ -s "$EGGNOG" ] || echo "WARNING: eggNOG file not found: $EGGNOG (will run without annotation)"

echo "Sample      : $SAMPLE"
echo "Isolate faa : $ISOLATE_FAA"
echo "References  : $REF_DIR"
echo ""

# --- gather reference proteomes, excluding the isolate itself ---
shopt -s nullglob globstar
REF_FAA=("$REF_DIR"/**/*.faa)
shopt -u nullglob globstar
[ "${#REF_FAA[@]}" -gt 0 ] || { echo "ERROR: no .faa under $REF_DIR" >&2; exit 1; }

: > refs_combined.faa
: > reference_genomes_used.txt
FIRST_ID=$(grep -m1 '^>' "$ISOLATE_FAA" | cut -d' ' -f1)
N_REF=0
for F in "${REF_FAA[@]}"; do
  # The isolate must never appear in its own reference set - every protein
  # would then have a perfect self-hit and nothing would look unique. Two
  # guards: same file, and same basename (catches the other Prokka run's copy).
  if [ "$(readlink -f "$F")" = "$(readlink -f "$ISOLATE_FAA")" ]; then
    echo "  excluded (same file)     : $F"; continue
  fi
  case "$(basename "$F")" in
    "$SAMPLE".faa) echo "  excluded (same basename) : $F"; continue ;;
  esac
  # Third guard: any path mentioning the sample name. Prokka output under a
  # different directory is often named for the assembly, not the sample, so
  # the two guards above can both miss the isolate's own copy.
  case "$F" in
    *"$SAMPLE"*) echo "  excluded (path contains $SAMPLE) : $F"; continue ;;
  esac
  # Fourth guard: identical first protein ID. Two Prokka runs of the SAME
  # genome share locus tags; two different genomes never do.
  if [ "$(grep -m1 '^>' "$F" | cut -d' ' -f1)" = "$FIRST_ID" ]; then
    echo "  excluded (same locus tags) : $F"; continue
  fi
  cat "$F" >> refs_combined.faa
  basename "$F" >> reference_genomes_used.txt
  N_REF=$((N_REF + 1))
done

N_QUERY=$(grep -c '^>' "$ISOLATE_FAA")
N_SUBJ=$(grep -c '^>' refs_combined.faa)
echo ""
echo "Isolate proteins   : $N_QUERY"
echo "Reference proteomes: $N_REF files, $N_SUBJ proteins"
[ "$N_REF" -ge 10 ] || echo "  !! only $N_REF references - report this n wherever you state 'unique'"
echo ""

# --- search ---
"$DIAMOND" makedb --in refs_combined.faa --db refs --threads "$SLURM_CPUS_PER_TASK" --quiet
"$DIAMOND" blastp \
  --query "$ISOLATE_FAA" \
  --db refs \
  --out hits.tsv \
  --outfmt 6 qseqid sseqid pident length evalue bitscore qcovhsp scovhsp \
  --evalue "$EVALUE" \
  --max-target-seqs 1 \
  --threads "$SLURM_CPUS_PER_TASK" \
  --sensitive \
  --quiet
echo "DIAMOND hits: $(wc -l < hits.tsv)"

# Self-inclusion detector. If the isolate slipped into its own reference set,
# nearly every protein gets a full-length 100% self-hit and the result is a
# silent zero. No real congener produces that distribution.
SELF=$(awk -F'\t' '$3 >= 99.9 && $7 >= 99 {n++} END{print n+0}' hits.tsv)
NQ=$(grep -c '^>' "$ISOLATE_FAA")
if [ "$NQ" -gt 0 ] && [ "$((100 * SELF / NQ))" -ge 80 ]; then
  echo "" >&2
  echo "ERROR: ${SELF}/${NQ} proteins have a full-length 100% identity hit." >&2
  echo "       The isolate is almost certainly inside its own reference set." >&2
  echo "       Check reference_genomes_used.txt and remove its own proteome." >&2
  exit 1
fi
echo "Full-length 100% hits: $SELF / $NQ  (high values mean self-inclusion)"

# --- classify and annotate ---
python3 - "$ISOLATE_FAA" hits.tsv "$EGGNOG" "$DIVERGENT_MAX_ID" "$SAMPLE" <<'PY'
import csv, sys, re
from collections import Counter

faa, hits_path, egg_path, div_max, name = sys.argv[1:6]
div_max = float(div_max)

queries = [l[1:].split()[0] for l in open(faa) if l.startswith(">")]

best = {}
for line in open(hits_path):
    p = line.rstrip("\n").split("\t")
    if len(p) < 8:
        continue
    q, pid, bits = p[0], float(p[2]), float(p[5])
    if q not in best or bits > best[q][1]:
        best[q] = (pid, bits, p[1])

egg, header = {}, None
try:
    for line in open(egg_path, errors="replace"):
        line = line.rstrip("\n")
        if line.startswith("##") or not line.strip():
            continue
        if line.startswith("#query"):
            header = line.lstrip("#").split("\t"); header[0] = "query"; continue
        if header is None:
            continue
        f = line.split("\t")
        if len(f) < len(header):
            f += [""] * (len(header) - len(f))
        d = dict(zip(header, f))
        egg[d["query"]] = d
except OSError:
    print("  !! eggNOG file unreadable; proceeding without functional annotation")

absent    = [q for q in queries if q not in best]
divergent = [q for q in queries if q in best and best[q][0] < div_max]
matched   = sum(1 for q in queries if q in egg)

print(f"\nProteins in isolate            : {len(queries)}")
print(f"  no homolog in any congener   : {len(absent)}  ({100*len(absent)/len(queries):.1f}%)")
print(f"  best hit below {div_max:.0f}% identity : {len(divergent)}  "
      f"({100*len(divergent)/len(queries):.1f}%)")
print(f"eggNOG identifier match        : {matched}/{len(queries)} "
      f"({100*matched/max(len(queries),1):.1f}%)")
if matched < 0.5 * len(queries):
    print("  !! LOW MATCH - ISOLATE_FAA is probably not the .faa you annotated.")
    print("     It must be the one under PROKKA_DIR/prokka_out, not the Roary")
    print("     Prokka run: identical gene calls, different locus tags.")

COG = {"J":"Translation, ribosomal structure and biogenesis","A":"RNA processing and modification",
 "K":"Transcription","L":"Replication, recombination and repair","B":"Chromatin structure",
 "D":"Cell cycle control, cell division","V":"Defense mechanisms","T":"Signal transduction",
 "M":"Cell wall/membrane/envelope biogenesis","N":"Cell motility","Z":"Cytoskeleton",
 "W":"Extracellular structures","U":"Intracellular trafficking, secretion",
 "O":"Posttranslational modification, protein turnover, chaperones",
 "X":"Mobilome: prophages, transposons","C":"Energy production and conversion",
 "G":"Carbohydrate transport and metabolism","E":"Amino acid transport and metabolism",
 "F":"Nucleotide transport and metabolism","H":"Coenzyme transport and metabolism",
 "I":"Lipid transport and metabolism","P":"Inorganic ion transport and metabolism",
 "Q":"Secondary metabolites biosynthesis, transport and catabolism",
 "R":"General function prediction only","S":"Function unknown","-":"Unassigned"}

def dump(tag, ids):
    rows, cog, hyp = [], Counter(), 0
    for q in ids:
        e = egg.get(q, {})
        c = (e.get("COG_category") or "-").strip() or "-"
        for letter in c:
            if letter in COG:
                cog[letter] += 1
        desc = e.get("Description", "")
        if not desc or re.search(r"hypothetical|unknown|uncharacter", desc, re.I):
            hyp += 1
        rows.append([q, e.get("Preferred_name",""), c,
                     "; ".join(COG.get(x,"") for x in c if x in COG),
                     desc, e.get("KEGG_ko",""), e.get("KEGG_Pathway",""),
                     f"{best[q][0]:.1f}" if q in best else "",
                     best[q][2] if q in best else ""])
    with open(f"{name}_{tag}.csv","w",newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Protein","Gene name","COG category","COG description",
                    "eggNOG description","KEGG KO","KEGG pathway",
                    "Best hit % identity","Best hit subject"])
        w.writerows(rows)
    tot = sum(v for k,v in cog.items() if k != "-")
    with open(f"{name}_{tag}_cog.csv","w",newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["COG category","Description","Proteins","% of COG-assigned"])
        for k,v in cog.most_common():
            w.writerow([k, COG[k], v, round(100*v/tot,1) if tot and k!="-" else ""])
    print(f"  wrote {name}_{tag}.csv and {name}_{tag}_cog.csv "
          f"({len(rows)} proteins, {hyp} unannotated/hypothetical = "
          f"{100*hyp/max(len(rows),1):.0f}%)")

print()
dump("absent_from_genus", absent)
dump("divergent", divergent)

# Best-hit identity for EVERY protein. The shape of this distribution shows how
# the isolate relates to its genus overall; its left tail is where the
# distinctive proteins are. Lets you re-cut the threshold without re-searching.
with open(f"{name}_best_hit_identity.csv","w",newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Protein","Best hit % identity","Best hit subject","Has homolog"])
    for q in queries:
        if q in best:
            w.writerow([q, f"{best[q][0]:.1f}", best[q][2], "yes"])
        else:
            w.writerow([q, "", "", "no"])
vals = sorted(best[q][0] for q in queries if q in best)
if vals:
    def pct(p):
        return vals[min(int(p / 100 * len(vals)), len(vals) - 1)]
    print(f"\n  Best-hit identity across {len(vals)} proteins with a homolog:")
    print(f"    median {pct(50):.1f}%   10th pct {pct(10):.1f}%   "
          f"1st pct {pct(1):.1f}%   min {vals[0]:.1f}%")
    for cut in (40, 50, 60, 70):
        n = sum(1 for v in vals if v < cut)
        print(f"    best hit < {cut}% identity: {n} proteins "
              f"({100*n/len(queries):.1f}% of proteome)")
    print(f"  wrote {name}_best_hit_identity.csv")
PY

echo ""
echo "Outputs in $OUTDIR:"
echo "  reference_genomes_used.txt        - the exact comparator set (n = $N_REF)"
echo "                                      cite this n in the text and list it"
echo "                                      as a supplementary table"
echo "  ${SAMPLE}_absent_from_genus.csv   - proteins with no homolog in any congener"
echo "  ${SAMPLE}_divergent.csv           - proteins whose best congener hit is weak"
echo "  ${SAMPLE}_best_hit_identity.csv   - best-hit identity for every protein"
echo "  *_cog.csv                         - COG profile of each set"
echo ""
echo "To run the other isolate: change SAMPLE at the top to Sample4 and resubmit."
echo ""
echo "CAVEAT: absence of a homolog in a DRAFT reference assembly can be an"
echo "assembly gap. Confirm any gene you highlight by tblastn against the"
echo "reference ASSEMBLIES (nucleotide), not just the proteomes."
