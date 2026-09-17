#!/bin/bash
#SBATCH --job-name=prokka_samples
#SBATCH -A YOUR_ACCOUNT_HERE
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -o prokka_samples_%j.out
#SBATCH -e prokka_samples_%j.err

set -euo pipefail

# Prokka annotation of your 6 sample genomes - the prerequisite for
# 13_funct_annotation.slurm (needs .faa), 14_antismash.slurm (prefers
# .gbk), and 15_mge_defense.slurm (DefenseFinder prefers .faa).
#
# Output layout matches exactly what those scripts expect:
#   $PROJECT/prokka_out/<Sample>/<Sample>.faa, <Sample>.gbk, <Sample>.gff, ...
#
# (--prefix is set to the sample name, so files are Sample1.faa etc., and
# the DefenseFinder input-matching in script 15 works without changes.)
#
# NOTE: this annotates ONLY the genomes sitting directly in $GENOME_DIR
# (your 6 samples) - NOT reference_genomes/. Your genus-wide reference
# Prokka runs for Roary are a separate thing and stay wherever you made
# them; scripts 13/14 deliberately only want the 6 samples (see README).

# ---- EDIT THESE ----
PROJECT=/project/YOUR_PROJECT_HERE
GENOME_DIR=$PROJECT            # your 6 sample FASTAs, directly here
# --------------------

OUTBASE=$PROJECT/prokka_out
mkdir -p "$OUTBASE"

# Prokka: try the Ceres module first (it's in the preinstalled list);
# fall back to a conda env at $PROJECT/envs/prokka_env if you made one.
if module load prokka 2>/dev/null && command -v prokka >/dev/null 2>&1; then
  PROKKA=prokka
  echo "Using Ceres prokka module: $(prokka --version 2>&1)"
elif [ -x "$PROJECT/envs/prokka_env/bin/prokka" ]; then
  export PATH="$PROJECT/envs/prokka_env/bin:$PATH"
  PROKKA=$PROJECT/envs/prokka_env/bin/prokka
  echo "Using conda prokka: $($PROKKA --version 2>&1)"
else
  echo "ERROR: prokka not found. Either 'module load prokka' doesn't exist" >&2
  echo "on this cluster (check: module spider prokka) or install it via" >&2
  echo "conda from a salloc session:" >&2
  echo "  mamba create -y --prefix $PROJECT/envs/prokka_env -c bioconda -c conda-forge prokka" >&2
  exit 1
fi

shopt -s nullglob
FA_FILES=("$GENOME_DIR"/*.fasta "$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fa)
shopt -u nullglob
if [ "${#FA_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no genome FASTAs found directly in $GENOME_DIR." >&2
  exit 1
fi
echo "Annotating ${#FA_FILES[@]} genomes."

for FA in "${FA_FILES[@]}"; do
  NAME=$(basename "$FA"); NAME=${NAME%.*}
  OUT=$OUTBASE/$NAME
  if [ -f "$OUT/$NAME.faa" ] && [ -f "$OUT/$NAME.gbk" ]; then
    echo "== $NAME already annotated, skipping =="
    continue
  fi
  echo "== Prokka: $NAME =="
  # --force: overwrite a half-finished dir from an interrupted run.
  # No --genus/--species set: with 6 different genera, per-sample genus
  # options would need a lookup table for marginal benefit - default
  # (kingdom Bacteria) annotation is what's normally used and reported.
  $PROKKA \
    --outdir "$OUT" \
    --prefix "$NAME" \
    --cpus "$SLURM_CPUS_PER_TASK" \
    --force \
    "$FA"
done

echo ""
echo "Done. Layout:"
find "$OUTBASE" -maxdepth 2 -name '*.faa' -o -maxdepth 2 -name '*.gbk' | sort
echo ""
echo "Now set PROKKA_DIR=$OUTBASE in scripts 13/14/15 (that's already their"
echo "default if PROJECT matches) and submit them:"
echo "  sbatch 13_funct_annotation.slurm"
echo "  sbatch 14_antismash.slurm"
echo "  sbatch 15_mge_defense.slurm"
