#!/bin/bash
#SBATCH --job-name=genus_placement
#SBATCH --output=genus_placement_%j.out
#SBATCH --error=genus_placement_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --partition=short
#SBATCH --account=kanglab

# ---------------------------------------------------------------------------
# Genus-level placement of Sample4: Massilia or Telluria?
#
# WHAT THIS DOES AND DOES NOT SETTLE
#
#   It does NOT decide the genus NAME. Massilia and Telluria are both validly
#   published genera. Bowman (2023, IJSEM 73:005991) reclassified 35 Massilia
#   species into Telluria; GTDB adopted that, which is why GTDB-Tk told you
#   "Telluria". LPSN did not adopt it: LPSN currently lists Massilia as the
#   correct name and Telluria haematophila as a homotypic SYNONYM of
#   Massilia haematophila. Choosing between them is a nomenclatural decision,
#   not a computation.
#
#   What this DOES settle is the empirical half, which you do need:
#     1. Is Sample4 congeneric with Massilia/Telluria at all?  (POCP > 50%)
#     2. Bowman split Massilia THREE ways - Telluria, Mokoshia and Zemynaea.
#        If you follow that scheme you must know which of the three Sample4
#        belongs to. GTDB said Telluria; this checks it independently.
#
# Metrics computed, against the type species of every candidate genus:
#     POCP  percentage of conserved proteins   (Qin et al. 2014) genus if > 50%
#     AAI   average amino-acid identity        genus if roughly > 65%
#     ANI   for reference only - it saturates below the genus level
#
# Submit:  sbatch 23_genus_placement.slurm
# ---------------------------------------------------------------------------

set -euo pipefail

PROJECT=/project/kanglab/rishi.bhandari
ENVS=$PROJECT/envs
OUTDIR=$PROJECT/genus_placement
REFDIR=$OUTDIR/type_genomes
THREADS=${SLURM_CPUS_PER_TASK:-16}

DIAMOND=$ENVS/funct_env/bin/diamond
SAMPLE=Sample4
SAMPLE_FAA=$PROJECT/PROKKA_DIR/prokka_out/$SAMPLE/$SAMPLE.faa

mkdir -p "$REFDIR" "$OUTDIR/faa" "$OUTDIR/work"
cd "$OUTDIR"

echo "=================================================================="
echo "Genus placement of $SAMPLE   $(date)"
echo "=================================================================="

if [ ! -s "$SAMPLE_FAA" ]; then
  echo "ERROR: proteome not found at $SAMPLE_FAA" >&2
  echo "This must be the SAME Prokka run used for the functional annotation." >&2
  exit 1
fi
echo "Isolate proteome: $SAMPLE_FAA ($(grep -c '^>' "$SAMPLE_FAA") proteins)"

# ---------------------------------------------------------------------------
# 1. Fetch the type-strain genomes that actually decide the question.
#
#    Massilia timonae      type species of Massilia
#    Telluria mixta        type species of Telluria
#    Massilia eurypsychrophila   type species of Mokoshia (Bowman 2023)
#    Massilia arenosa            the sole Zemynaea species (Bowman 2023)
#    Massilia haematophila       your closest relative by dDDH (51.3%)
# ---------------------------------------------------------------------------
echo
echo "== Step 1: downloading type-strain genomes =="

# Prefer the datasets CLI if it is on PATH; adjust if yours lives elsewhere.
DATASETS=$(command -v datasets || echo "$ENVS/funct_env/bin/datasets")
if [ ! -x "$DATASETS" ]; then
  echo "ERROR: NCBI 'datasets' CLI not found. Load it or install into an env." >&2
  exit 1
fi

TAXA=(
  "Massilia timonae"
  "Telluria mixta"
  "Massilia eurypsychrophila"
  "Massilia arenosa"
  "Massilia haematophila"
)

for T in "${TAXA[@]}"; do
  SLUG=$(echo "$T" | tr ' ' '_')
  if [ -s "$REFDIR/$SLUG/protein.faa" ]; then
    echo "  $T already present, skipping"
    continue
  fi
  echo "  fetching $T"
  mkdir -p "$REFDIR/$SLUG"
  # --reference is deliberately NOT used: type strains are often not the
  # designated reference assembly. Take the type-material assembly instead.
  "$DATASETS" download genome taxon "$T" \
      --assembly-source RefSeq \
      --include genome,protein \
      --filename "$REFDIR/$SLUG.zip" || {
        echo "  !! download failed for $T - continuing" >&2; continue; }
  unzip -qo "$REFDIR/$SLUG.zip" -d "$REFDIR/$SLUG/raw"
  # Take the first assembly that carries a protein set.
  FAA=$(find "$REFDIR/$SLUG/raw" -name "protein.faa" | head -n1)
  FNA=$(find "$REFDIR/$SLUG/raw" -name "*_genomic.fna" | head -n1)
  if [ -n "$FAA" ]; then cp "$FAA" "$REFDIR/$SLUG/protein.faa"; fi
  if [ -n "$FNA" ]; then cp "$FNA" "$REFDIR/$SLUG/genomic.fna"; fi
done

echo
echo "  IMPORTANT - verify type material by hand before trusting the result:"
echo "    datasets summary genome taxon 'Massilia timonae' --as-json-lines \\"
echo "      | dataformat tsv genome --fields accession,organism-name,assminfo-biosample-strain,checkm-completeness"
echo "  and confirm the accession you used is the one NCBI flags as"
echo "  'assembly from type material'. Getting a non-type strain here would"
echo "  invalidate the whole comparison."

# ---------------------------------------------------------------------------
# 2. Assemble the proteome set
# ---------------------------------------------------------------------------
echo
echo "== Step 2: collecting proteomes =="
cp "$SAMPLE_FAA" "$OUTDIR/faa/$SAMPLE.faa"
for T in "${TAXA[@]}"; do
  SLUG=$(echo "$T" | tr ' ' '_')
  if [ -s "$REFDIR/$SLUG/protein.faa" ]; then
    cp "$REFDIR/$SLUG/protein.faa" "$OUTDIR/faa/$SLUG.faa"
    echo "  $SLUG: $(grep -c '^>' "$OUTDIR/faa/$SLUG.faa") proteins"
  else
    echo "  !! $SLUG has no proteome - it will be missing from the table" >&2
  fi
done

# ---------------------------------------------------------------------------
# 3. All-vs-all DIAMOND, both directions (POCP and AAI are asymmetric inputs)
# ---------------------------------------------------------------------------
echo
echo "== Step 3: reciprocal DIAMOND searches =="
cd "$OUTDIR/work"

for F in "$OUTDIR"/faa/*.faa; do
  B=$(basename "$F" .faa)
  [ -s "$B.dmnd" ] || "$DIAMOND" makedb --in "$F" --db "$B" --quiet
done

for Q in "$OUTDIR"/faa/*.faa; do
  QB=$(basename "$Q" .faa)
  for S in "$OUTDIR"/faa/*.faa; do
    SB=$(basename "$S" .faa)
    [ "$QB" = "$SB" ] && continue
    OUT="${QB}__vs__${SB}.tsv"
    [ -s "$OUT" ] && continue
    echo "  $QB -> $SB"
    "$DIAMOND" blastp \
      --query "$Q" --db "$SB" --out "$OUT" \
      --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
      --evalue 1e-5 --max-target-seqs 1 --threads "$THREADS" \
      --sensitive --quiet
  done
done

# ---------------------------------------------------------------------------
# 4. Compute POCP and AAI
# ---------------------------------------------------------------------------
echo
echo "== Step 4: POCP and AAI =="
python3 "$PROJECT/scripts/24_pocp_aai.py" \
    --faa-dir "$OUTDIR/faa" \
    --hits-dir "$OUTDIR/work" \
    --focal "$SAMPLE" \
    --out "$OUTDIR/genus_placement_results.csv"

echo
echo "=================================================================="
echo "Done. Results: $OUTDIR/genus_placement_results.csv"
echo
echo "How to read it:"
echo "  POCP > 50%  => same genus as that type species (Qin et al. 2014)"
echo "  AAI  > ~65% => same genus, corroborating"
echo "  Compare Sample4's POCP against Massilia timonae vs Telluria mixta"
echo "  vs Massilia eurypsychrophila. The highest, and any above 50%, tell"
echo "  you which genus concept Sample4 sits in."
echo
echo "Then, separately, decide which NAME to publish under. See the note at"
echo "the top of this script - that part is a choice, not a measurement."
echo "=================================================================="
