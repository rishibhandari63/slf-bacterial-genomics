#!/bin/bash
#SBATCH --job-name=gtdbtk_classify
#SBATCH -A YOUR_ACCOUNT_HERE
#SBATCH --partition=mem
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH -o gtdbtk_classify_%j.out
#SBATCH -e gtdbtk_classify_%j.err

set -euo pipefail

# ---- EDIT THESE ----
PROJECT=/project/YOUR_PROJECT_HERE
WORKDIR=$PROJECT/genome_phylogeny
GENOME_DIR=$WORKDIR/genomes
GENOME_EXT=fna          # change to fasta/fa if that's what your files use
# ---------------------

OUTDIR=$WORKDIR/gtdbtk_classify_out
mkdir -p "$OUTDIR"

module load gtdbtk

# Point at the database you downloaded in step 1, unless the module already
# ships with GTDBTK_DATA_PATH preset (check with: echo $GTDBTK_DATA_PATH
# right after `module load gtdbtk`, before this override).
if [ -f "$PROJECT/gtdbtk_data/env.sh" ]; then
  source "$PROJECT/gtdbtk_data/env.sh"
fi
echo "Using GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH:-<not set - will fail>}"

# Note: no --mash_db (or --skani_db) flag here. Current gtdbtk (2.4.x on
# Ceres) does ANI pre-screening with skani, not Mash, and handles that
# internally against GTDBTK_DATA_PATH - there is no separate cache file
# to point it at. Passing --mash_db on this version errors out with
# "unrecognized arguments".
gtdbtk classify_wf \
  --genome_dir "$GENOME_DIR" \
  --extension "$GENOME_EXT" \
  --out_dir "$OUTDIR" \
  --cpus "$SLURM_CPUS_PER_TASK" \
  --pplacer_cpus 4

echo ""
echo "Key outputs:"
echo "  Taxonomy + ANI to closest reference: $OUTDIR/gtdbtk.bac120.summary.tsv"
echo "  Placement tree (Newick):             $OUTDIR/gtdbtk.bac120.classify.tree"
echo "  (If any genome is Archaea instead, look for ar53.* equivalents)"
echo ""
echo "Open the summary.tsv and note, for each genome: classification,"
echo "closest_genome_reference, and closest_genome_ani. An ANI >= ~95%"
echo "against a named species usually means that's your species; lower"
echo "values mean it's likely a novel species within that genus, and the"
echo "classify.tree tells you which established lineage it's nested in."
