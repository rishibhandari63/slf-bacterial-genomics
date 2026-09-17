#!/bin/bash
#SBATCH --job-name=funct_annot
#SBATCH -A YOUR_ACCOUNT_HERE
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o funct_annot_%j.out
#SBATCH -e funct_annot_%j.err

set -euo pipefail

# Functional annotation of your 6 isolates: eggNOG-mapper (COG categories,
# KEGG KOs/pathways, gene names) + run_dbcan v5 (CAZymes).
#
# INPUT: protein FASTAs (.faa) from your Prokka runs - the same
# annotations you fed Roary. Point PROKKA_DIR at the parent folder that
# holds one Prokka output dir per sample; the glob below finds each
# sample's .faa wherever Prokka put it.

# ---- EDIT THESE ----
PROJECT=/project/YOUR_PROJECT_HERE
# Parent dir containing your Prokka output folders for the 6 samples
# (e.g. $PROKKA_DIR/Sample1/Sample1.faa or similar layout):
PROKKA_DIR=$PROJECT/prokka_out
# --------------------

ENVS=$PROJECT/envs
DBS=$PROJECT/functional_dbs
OUTDIR=$PROJECT/funct_annotation_out
mkdir -p "$OUTDIR"

# eggNOG-mapper needs its env's bin on PATH for its internal diamond calls
# (same subprocess-PATH lesson we learned the hard way with pyani/nucmer).
export PATH="$ENVS/funct_env/bin:$PATH"

shopt -s nullglob
FAA_FILES=("$PROKKA_DIR"/*/*.faa "$PROKKA_DIR"/*.faa)
shopt -u nullglob
if [ "${#FAA_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no .faa files found under $PROKKA_DIR (or one level down)." >&2
  echo "Point PROKKA_DIR at the folder holding your Prokka outputs." >&2
  exit 1
fi
echo "Found ${#FAA_FILES[@]} protein FASTAs:"
printf '  %s\n' "${FAA_FILES[@]}"

for FAA in "${FAA_FILES[@]}"; do
  NAME=$(basename "$FAA"); NAME=${NAME%.faa}

  # --- eggNOG-mapper ---
  EGG_OUT=$OUTDIR/eggnog/$NAME
  if [ -f "$EGG_OUT/${NAME}.emapper.annotations" ]; then
    echo "== eggNOG: $NAME already done, skipping =="
  else
    echo "== eggNOG: $NAME =="
    mkdir -p "$EGG_OUT"
    "$ENVS/funct_env/bin/emapper.py" \
      -i "$FAA" \
      --itype proteins \
      -o "$NAME" \
      --output_dir "$EGG_OUT" \
      --data_dir "$DBS/eggnog" \
      --cpu "$SLURM_CPUS_PER_TASK" \
      --override
  fi

  # --- dbCAN (CAZymes) ---
  CAZY_OUT=$OUTDIR/dbcan/$NAME
  if [ -d "$CAZY_OUT" ] && [ -n "$(ls -A "$CAZY_OUT" 2>/dev/null)" ]; then
    echo "== dbCAN: $NAME already done, skipping =="
  else
    echo "== dbCAN: $NAME =="
    mkdir -p "$CAZY_OUT"
    "$ENVS/funct_env/bin/run_dbcan" CAZyme_annotation \
      --input_raw_data "$FAA" \
      --mode protein \
      --output_dir "$CAZY_OUT" \
      --db_dir "$DBS/dbcan" \
      --threads "$SLURM_CPUS_PER_TASK"
  fi
done

echo ""
echo "Done. Key outputs per sample:"
echo "  $OUTDIR/eggnog/<sample>/<sample>.emapper.annotations"
echo "    - one row per protein: COG category, KEGG KO/pathway, gene name, description."
echo "  $OUTDIR/dbcan/<sample>/  - CAZyme family calls (overview table = consensus of HMM/DIAMOND/dbCAN-sub;"
echo "    standard practice: count a protein as a CAZyme when >=2 of 3 tools agree)."
echo ""
echo "For the comparative angle: join eggnog annotations against Roary's"
echo "gene_presence_absence.csv to functionally profile genes unique to your"
echo "isolates vs. their genus relatives."
