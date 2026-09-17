#!/bin/bash
#SBATCH --job-name=checkm2_qc
#SBATCH -A YOUR_ACCOUNT_HERE
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o checkm2_%j.out
#SBATCH -e checkm2_%j.err

set -euo pipefail

# ---- EDIT THESE ----
PROJECT=/project/YOUR_PROJECT_HERE
WORKDIR=$PROJECT/genome_phylogeny
GENOME_DIR=$WORKDIR/genomes
GENOME_EXT=fna          # change to fasta/fa if that's what your files use
# ---------------------

CHECKM2_DB_DIR=$PROJECT/checkm2_db
OUTDIR=$WORKDIR/checkm2_out
mkdir -p "$OUTDIR"

module load checkm2

# Download the CheckM2 DIAMOND reference DB once, reuse afterwards.
# Check specifically for the .dmnd file rather than "dir non-empty" - a
# previous run can crash after extraction but leave a populated directory.
#
# --no_write_json_db: the checkm2 module lives in a shared, read-only
# /software install, so checkm2 database --download will crash with a
# PermissionError on its last step (writing a bookkeeping JSON pointer
# file back into its own install dir) unless this flag is set. The actual
# database download/extraction/checksum happens BEFORE that step, so if
# you already hit that crash once, check for an existing .dmnd file below
# before assuming you need to redownload ~all of it again.
if [ -z "$(find "$CHECKM2_DB_DIR" -name '*.dmnd' 2>/dev/null)" ]; then
  echo "CheckM2 database not found at $CHECKM2_DB_DIR, downloading..."
  mkdir -p "$CHECKM2_DB_DIR"
  checkm2 database --download --path "$CHECKM2_DB_DIR" --no_write_json_db
else
  echo "Found existing CheckM2 database under $CHECKM2_DB_DIR, skipping download."
fi

checkm2 predict \
  --threads "$SLURM_CPUS_PER_TASK" \
  --input "$GENOME_DIR" \
  --extension "$GENOME_EXT" \
  --output-directory "$OUTDIR" \
  --database_path "$(find "$CHECKM2_DB_DIR" -name '*.dmnd' | head -n1)"

echo "CheckM2 results: $OUTDIR/quality_report.tsv"
echo "Check completeness (>~90% is typically fine) and contamination"
echo "(<~5% is typically fine) before trusting the phylogeny downstream."
