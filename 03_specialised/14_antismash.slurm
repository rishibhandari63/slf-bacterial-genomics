#!/bin/bash
#SBATCH --job-name=antismash
#SBATCH -A YOUR_ACCOUNT_HERE
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o antismash_%j.out
#SBATCH -e antismash_%j.err

set -euo pipefail

# Secondary metabolite / biosynthetic gene cluster (BGC) detection with
# antiSMASH for your 6 isolates.
#
# INPUT: Prokka's .gbk output per sample if available (preferred - keeps
# your existing gene calls and locus tags, so BGC genes are directly
# cross-referenceable with the eggNOG/Roary tables), otherwise falls back
# to the raw genome FASTA with antiSMASH's own gene calling (prodigal).

# ---- EDIT THESE ----
PROJECT=/project/YOUR_PROJECT_HERE
PROKKA_DIR=$PROJECT/prokka_out          # same layout as script 13
GENOME_DIR=$PROJECT                     # raw sample FASTAs (fallback input)
# --------------------

ENVS=$PROJECT/envs
DBS=$PROJECT/functional_dbs
OUTDIR=$PROJECT/antismash_out
mkdir -p "$OUTDIR"

export PATH="$ENVS/antismash_env/bin:$PATH"
echo "antismash: $("$ENVS/antismash_env/bin/antismash" --version)"

# Prefer Prokka .gbk; fall back to raw FASTA per-sample.
shopt -s nullglob
GBK_FILES=("$PROKKA_DIR"/*/*.gbk "$PROKKA_DIR"/*.gbk)
shopt -u nullglob

declare -A INPUTS=()
if [ "${#GBK_FILES[@]}" -gt 0 ]; then
  echo "Using ${#GBK_FILES[@]} Prokka GenBank files as input."
  for GBK in "${GBK_FILES[@]}"; do
    NAME=$(basename "$GBK"); NAME=${NAME%.gbk}
    INPUTS[$NAME]=$GBK
  done
else
  echo "No Prokka .gbk files found under $PROKKA_DIR - falling back to raw"
  echo "genome FASTAs in $GENOME_DIR (antiSMASH will run prodigal itself)."
  shopt -s nullglob
  FA_FILES=("$GENOME_DIR"/*.fasta "$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fa)
  shopt -u nullglob
  if [ "${#FA_FILES[@]}" -eq 0 ]; then
    echo "ERROR: no inputs found in $PROKKA_DIR or $GENOME_DIR." >&2
    exit 1
  fi
  for FA in "${FA_FILES[@]}"; do
    NAME=$(basename "$FA"); NAME=${NAME%.*}
    INPUTS[$NAME]=$FA
  done
fi

for NAME in "${!INPUTS[@]}"; do
  SRC=${INPUTS[$NAME]}
  SAMPLE_OUT=$OUTDIR/$NAME
  if [ -f "$SAMPLE_OUT/index.html" ]; then
    echo "== $NAME already done, skipping =="
    continue
  fi
  rm -rf "$SAMPLE_OUT"
  echo "== antiSMASH: $NAME ($SRC) =="
  EXTRA_ARGS=()
  case "$SRC" in
    *.gbk) : ;;  # annotated input - use existing gene calls
    *) EXTRA_ARGS+=(--genefinding-tool prodigal) ;;
  esac
  "$ENVS/antismash_env/bin/antismash" \
    --databases "$DBS/antismash" \
    --output-dir "$SAMPLE_OUT" \
    --cpus "$SLURM_CPUS_PER_TASK" \
    --cb-general --cb-knownclusters --cb-subclusters \
    --pfam2go \
    "${EXTRA_ARGS[@]}" \
    "$SRC"
done

echo ""
echo "Done. Per sample: $OUTDIR/<sample>/index.html is the interactive"
echo "report (download the whole sample folder and open locally in a"
echo "browser). Region GenBank files (*.region*.gbk) are the per-BGC"
echo "records. --cb-knownclusters output tells you % similarity of each"
echo "cluster to characterized MIBiG clusters - clusters with low/no"
echo "similarity to anything known are the interesting ones for a"
echo "'novel BGC' claim."
