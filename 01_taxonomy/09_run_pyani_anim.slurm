#!/bin/bash
#SBATCH --job-name=pyani_anim
#SBATCH -A YOUR_ACCOUNT_HERE
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -o pyani_anim_%j.out
#SBATCH -e pyani_anim_%j.err

set -euo pipefail

# ---- EDIT THESE ----
PROJECT=/project/YOUR_PROJECT_HERE
WORKDIR=$PROJECT                       # everything lives directly under $PROJECT
GENOME_DIR=$WORKDIR                    # your 6 query genomes (Sample1.fasta etc, directly here)
REF_DIR=$WORKDIR/reference_genomes     # populated by 05_download_refs_and_fastani.slurm
# No GENOME_EXT needed - both loops below match .fasta/.fna/.fa automatically
# and standardize everything to .fna inside the scratch directory, so this
# script no longer cares what extension your source files actually use.

# Display label for every genome that will appear in the figure. Key = the
# genome file's basename WITHOUT extension (this must match exactly - it's
# how pyani's --labels/--classes files join back to each genome). Covers
# all ~34 genomes in reference_genomes/ (every reference is now included,
# not just the closest 6 - see the glob below) plus your 6 samples; any
# genome not listed here just falls back to displaying its raw accession.
declare -A LABEL_MAP=(
  [Sample1]="Acinetobacter schindleri (Sample1)"
  [Sample2]="Chryseobacterium takakiae (Sample2)"
  [Sample3]="Sphingomonas yabuuchiae (Sample3)"
  [Sample4]="Telluria sp. (Sample4)"
  [Sample5]="Glutamicibacter mishrai (Sample5)"
  [Sample8]="Agromyces sp. (Sample8)"

  [GCF_000368625.1]="A. schindleri (ref, ANI 97.2)"
  [GCF_022700675.1]="Acinetobacter sp. (ANI 91.8)"
  [GCF_002688565.1]="Acinetobacter sp. (ANI 89.2)"
  [GCA_000761495.1]="A. idrijaensis (ANI 85.9)"
  [GCF_019343495.1]="A. lwoffii_P (ANI 85.6)"
  [GCF_000369625.1]="A. variabilis (ANI 84.9)"

  [GCF_900129385.1]="C. takakiae (ref, ANI 99.0)"
  [GCF_023195755.1]="C. nepalense (ANI 90.2)"
  [GCA_913778275.1]="Chryseobacterium sp. (ANI 88.6)"
  [GCF_003815835.1]="C. shandongense (ANI 87.6)"
  [GCF_900182825.1]="C. profundimaris (ANI 87.4)"
  [GCF_000708615.2]="C. hispalense (ANI 87.1)"

  [GCF_014199595.1]="S. yabuuchiae (ref, ANI 96.8)"
  [GCF_001477495.1]="S. yabuuchiae_A (ANI 92.9)"
  [GCF_000787715.1]="S. parapaucimobilis (ANI 91.9)"
  [GCF_019297835.1]="S. sanguinis_D (ANI 91.8)"
  [GCA_019112965.1]="S. excrementigallinarum (ANI 90.6)"
  [GCF_001591005.1]="S. sanguinis (ANI 89.9)"

  [GCA_003484545.1]="Telluria sp003484545 (ref, ANI 96.6)"
  [GCA_003520365.1]="Telluria sp. (ANI 94.5)"
  [GCA_002354135.1]="Telluria sp. (ANI 91.7)"
  [GCF_030411995.1]="Telluria sp. (ANI 87.2)"
  [GCF_020708045.1]="Telluria sp. (ANI 87.1)"
  [GCF_030711405.1]="T. brevitalea (ANI 86.8)"

  [GCF_012221945.1]="G. mishrai (ref, ANI 98.7)"
  [GCF_025854375.1]="Glutamicibacter sp. (ANI 92.2)"
  [GCF_900094805.1]="Glutamicibacter sp. (ANI 87.0)"
  [GCF_001302565.1]="G. halophytocola (ANI 86.7)"
  [GCF_006539525.1]="G. nicotianae (ANI 84.2)"

  [GCF_900113185.1]="Agromyces sp900113185 (ref, ANI 91.7)"
  [GCF_009749385.1]="A. allii (ANI 89.3)"
  [GCF_001660485.1]="A. aureus (ANI 88.2)"
  [GCF_001421565.1]="Agromyces sp. (ANI 87.5)"
  [GCF_009749465.1]="A. terreus (ANI 85.8)"
  [GCF_009749395.1]="A. salentinus (ANI 84.9)"
)

# Genus (or any grouping you want as the colored side-bar in the seaborn
# clustermap) - same keys as LABEL_MAP. Anything not listed falls back to
# "unknown" rather than erroring.
declare -A CLASS_MAP=(
  [Sample1]="Acinetobacter" [Sample2]="Chryseobacterium" [Sample3]="Sphingomonas"
  [Sample4]="Telluria" [Sample5]="Glutamicibacter" [Sample8]="Agromyces"

  [GCF_000368625.1]="Acinetobacter" [GCF_022700675.1]="Acinetobacter"
  [GCF_002688565.1]="Acinetobacter" [GCA_000761495.1]="Acinetobacter"
  [GCF_019343495.1]="Acinetobacter" [GCF_000369625.1]="Acinetobacter"

  [GCF_900129385.1]="Chryseobacterium" [GCF_023195755.1]="Chryseobacterium"
  [GCA_913778275.1]="Chryseobacterium" [GCF_003815835.1]="Chryseobacterium"
  [GCF_900182825.1]="Chryseobacterium" [GCF_000708615.2]="Chryseobacterium"

  [GCF_014199595.1]="Sphingomonas" [GCF_001477495.1]="Sphingomonas"
  [GCF_000787715.1]="Sphingomonas" [GCF_019297835.1]="Sphingomonas"
  [GCA_019112965.1]="Sphingomonas" [GCF_001591005.1]="Sphingomonas"

  [GCA_003484545.1]="Telluria" [GCA_003520365.1]="Telluria"
  [GCA_002354135.1]="Telluria" [GCF_030411995.1]="Telluria"
  [GCF_020708045.1]="Telluria" [GCF_030711405.1]="Telluria"

  [GCF_012221945.1]="Glutamicibacter" [GCF_025854375.1]="Glutamicibacter"
  [GCF_900094805.1]="Glutamicibacter" [GCF_001302565.1]="Glutamicibacter"
  [GCF_006539525.1]="Glutamicibacter"

  [GCF_900113185.1]="Agromyces" [GCF_009749385.1]="Agromyces"
  [GCF_001660485.1]="Agromyces" [GCF_001421565.1]="Agromyces"
  [GCF_009749465.1]="Agromyces" [GCF_009749395.1]="Agromyces"
)
# ---------------------

ENV_PREFIX=$PROJECT/envs/pyani_env
# pyani calls `average_nucleotide_identity.py` fine when invoked by full
# path (below), but INTERNALLY it spawns nucmer/delta-filter/show-coords
# as subprocesses that resolve by bare command name via $PATH - and since
# we deliberately skip `source activate`/`conda activate` (broken on this
# cluster), $PATH never picked up this env's bin/, where the pinned
# mummer=3.23 actually lives. Without this, those subprocess calls either
# fail to find nucmer/delta-filter at all, or fall through to some other
# MUMmer on the default PATH - which is exactly what caused "NUCmer
# comparison failed" / "No filter files found" in the last run. Exporting
# PATH directly (rather than via any activation mechanism) is the
# reliable fix regardless of what's wrong with activation here.
export PATH="$ENV_PREFIX/bin:$PATH"
echo "nucmer resolved to: $(which nucmer)"
"$ENV_PREFIX/bin/nucmer" --version
INDIR=$WORKDIR/pyani_input
OUTDIR=$WORKDIR/pyani_anim_out
rm -rf "$INDIR"
mkdir -p "$INDIR"
# Deliberately NOT pre-creating $OUTDIR: average_nucleotide_identity.py
# insists on creating its own -o directory and errors out ("exists
# (exiting)") if it already exists, even empty - that's its own safeguard
# against silently overwriting a previous run. -f/--force below also lets
# you rerun this script over a previous attempt without deleting anything
# by hand.
rm -rf "$OUTDIR"

# pyani wants every input genome as its own FASTA in one flat directory.
# Symlinking (not copying) keeps this cheap and re-runnable. Extension-
# agnostic and always renamed to .fna in $INDIR, same as the reference
# loop below - your query genomes are directly in $GENOME_DIR (not a
# genomes/ subfolder) and have been renamed .fna -> .fasta at least once
# already, so don't assume any particular extension here.
shopt -s nullglob
QUERY_FILES=("$GENOME_DIR"/*.fasta "$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fa)
shopt -u nullglob
if [ "${#QUERY_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no genomes found directly in $GENOME_DIR - check the path." >&2
  exit 1
fi
for SRC in "${QUERY_FILES[@]}"; do
  STEM=$(basename "$SRC")
  STEM=${STEM%.*}
  ln -sf "$SRC" "$INDIR/${STEM}.fna"
done
echo "Linked ${#QUERY_FILES[@]} query genomes from $GENOME_DIR."

# Every genome in $REF_DIR goes in - not just a curated subset. Same
# extension-agnostic glob + standardize-to-.fna pattern as the query loop.
shopt -s nullglob
REF_FILES=("$REF_DIR"/*.fasta "$REF_DIR"/*.fna "$REF_DIR"/*.fa)
shopt -u nullglob
if [ "${#REF_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no genomes found in $REF_DIR - check the path." >&2
  exit 1
fi
for SRC in "${REF_FILES[@]}"; do
  STEM=$(basename "$SRC")
  STEM=${STEM%.*}
  ln -sf "$SRC" "$INDIR/${STEM}.fna"
done
echo "Linked ${#REF_FILES[@]} reference genomes from $REF_DIR."

echo "Genomes going into ANIm ($(find "$INDIR" -name '*.fna' | wc -l) total):"
find "$INDIR" -name '*.fna' -exec basename {} \;

# Build pyani's --labels / --classes files. Format (confirmed from pyani's
# own source, pyani_tools.get_labels): three tab-separated columns,
# <anything>\t<genome file stem>\t<display text> - the first column is
# read and discarded, so a placeholder is fine; the SECOND column is the
# actual join key and must exactly equal the genome's filename without
# extension.
LABELS_FILE="$INDIR/labels.tab"
CLASSES_FILE="$INDIR/classes.tab"
: > "$LABELS_FILE"
: > "$CLASSES_FILE"
for f in "$INDIR"/*.fna; do
  STEM=$(basename "$f" ".fna")
  LABEL="${LABEL_MAP[$STEM]:-$STEM}"
  CLASS="${CLASS_MAP[$STEM]:-unknown}"
  printf 'placeholder\t%s\t%s\n' "$STEM" "$LABEL" >> "$LABELS_FILE"
  printf 'placeholder\t%s\t%s\n' "$STEM" "$CLASS" >> "$CLASSES_FILE"
done

# Call the pyani binary directly by full path rather than activating the
# env - `source activate` proved unreliable on this cluster earlier in
# this project, so we skip it entirely here.
# Log file deliberately points OUTSIDE $OUTDIR (into $INDIR, which already
# exists) - $OUTDIR must not exist yet when pyani starts, so a log path
# inside it would fail to open before pyani gets a chance to create it.
"$ENV_PREFIX/bin/average_nucleotide_identity.py" \
  -i "$INDIR" \
  -o "$OUTDIR" \
  -m ANIm \
  -g --gformat png,pdf,svg --gmethod seaborn \
  --labels "$LABELS_FILE" \
  --classes "$CLASSES_FILE" \
  --workers "$SLURM_CPUS_PER_TASK" \
  -v -l "$INDIR/pyani_run.log"

echo ""
echo "Done. Look for:"
echo "  $OUTDIR/ANIm_percentage_identity.png   <- the figure you're after"
echo "  $OUTDIR/ANIm_percentage_identity.tab   <- the underlying matrix"
echo "  $OUTDIR/ANIm_alignment_coverage.png    <- how much of each genome pair actually aligned"
echo ""
echo "scp ANIm_percentage_identity.png back to your laptop to view it."
echo ""
echo "Note: if nucmer errors out, double check 08_install_pyani_env.sh"
echo "actually installed mummer=3.23 into the env (run"
echo "  $ENV_PREFIX/bin/nucmer --version"
echo "and confirm it reports 3.23, not the Ceres 'mummer' module's 4.x)."
