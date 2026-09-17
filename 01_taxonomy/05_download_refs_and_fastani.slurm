#!/bin/bash
#SBATCH --job-name=fastani
#SBATCH -A YOUR_ACCOUNT_HERE
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -o fastani_%j.out
#SBATCH -e fastani_%j.err

set -euo pipefail

# ---- EDIT THESE ----
PROJECT=/project/YOUR_PROJECT_HERE
WORKDIR=$PROJECT/genome_phylogeny
GENOME_DIR=$WORKDIR/genomes
GENOME_EXT=fna

# Fill this in AFTER looking at step 3's gtdbtk.bac120.summary.tsv:
# one NCBI GenBank/RefSeq assembly accession per line, for the closest
# relative(s) of each of your genomes (and/or relevant type strains).
# Example accessions shown are placeholders - replace them.
REF_ACCESSIONS=(
  GCF_000000000.1
  GCF_000000000.2
)
# ---------------------

REF_DIR=$WORKDIR/reference_genomes
OUTDIR=$WORKDIR/fastani_out
mkdir -p "$REF_DIR" "$OUTDIR"

ENV_PREFIX=$PROJECT/envs/fastani_env
module load miniconda
source activate "$ENV_PREFIX"

# Download reference genomes via NCBI datasets CLI if not already present.
for ACC in "${REF_ACCESSIONS[@]}"; do
  TARGET="$REF_DIR/${ACC}.fna"
  if [ -f "$TARGET" ]; then
    echo "$ACC already downloaded, skipping."
    continue
  fi
  echo "Downloading $ACC..."
  datasets download genome accession "$ACC" \
    --include genome \
    --filename "$REF_DIR/${ACC}.zip"
  unzip -o -q "$REF_DIR/${ACC}.zip" -d "$REF_DIR/${ACC}_extract"
  # NCBI datasets nests the fasta a few directories deep; find and flatten it.
  find "$REF_DIR/${ACC}_extract" -name "*.fna" -exec cp {} "$TARGET" \;
  rm -rf "$REF_DIR/${ACC}.zip" "$REF_DIR/${ACC}_extract"
done

# Build query and reference lists for all-vs-all fastANI.
find "$GENOME_DIR" -name "*.${GENOME_EXT}" > "$OUTDIR/query_list.txt"
find "$REF_DIR" -name "*.fna" > "$OUTDIR/ref_list.txt"

fastANI \
  --ql "$OUTDIR/query_list.txt" \
  --rl "$OUTDIR/ref_list.txt" \
  -o "$OUTDIR/fastani_results.tsv" \
  -t "$SLURM_CPUS_PER_TASK"

echo ""
echo "Results: $OUTDIR/fastani_results.tsv"
echo "Columns: query, reference, ANI (%), count of orthologous fragments, total fragments."
echo "~95-96% ANI against a named type-strain genome = likely the same species."
echo "For a formal species-novelty claim, also consider running your"
echo "genomes through TYGS (https://tygs.dsmz.de/) which adds digital DDH"
echo "and a properly formatted phylogenomic tree for a species description."
