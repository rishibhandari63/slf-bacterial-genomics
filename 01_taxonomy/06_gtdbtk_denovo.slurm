#!/bin/bash
#SBATCH --job-name=gtdbtk_denovo
#SBATCH -A YOUR_ACCOUNT_HERE
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH -o gtdbtk_denovo_%j.out
#SBATCH -e gtdbtk_denovo_%j.err

set -euo pipefail

# OPTIONAL step. classify_wf (script 03) already gives you a placement
# tree; run this if you want an actual de novo ML tree built from your
# genomes + a small, hand-picked set of relatives - NOT the full GTDB
# reference tree.
#
# IMPORTANT, checked directly against the current GTDB-Tk docs: --taxa_filter
# (used in an earlier version of this script) only supports filtering at the
# DOMAIN or PHYLUM rank - there is no family/order-level filter. Your 6
# genomes span 3 phyla (Pseudomonadota, Bacteroidota, Actinomycetota), each
# of which has tens of thousands of reference genomes in GTDB - a
# phylum-filtered de_novo_wf run is not a SLURM-job-sized task, it would
# need distributed-computing-scale resources.
#
# So this script uses the documented workaround instead: --skip_gtdb_refs
# builds the tree from ONLY the genomes you hand it - your 6 queries plus
# every reference genome already sitting in $REF_DIR (the same folder used
# for the pyani/fastANI steps - whatever's in there goes in, no need to
# list accessions by hand) - plus one outgroup genome this script
# downloads automatically, picked from the GTDB taxonomy file you already
# have from step 1. With ~34 references that's ~41 genomes total - still a
# completely normal single-node job, just a bit slower than the earlier
# 13-genome version.

# ---- EDIT THESE ----
PROJECT=/project/YOUR_PROJECT_HERE
WORKDIR=$PROJECT                        # everything lives directly under $PROJECT
GENOME_DIR=$WORKDIR                     # your 6 query genomes (Sample1.fasta etc, directly here)
REF_DIR=$WORKDIR/reference_genomes      # from 05_download_refs_and_fastani.slurm
# No GENOME_EXT needed - genome discovery below matches .fasta/.fna/.fa
# automatically and standardizes everything to .fna in the scratch
# directory, so gtdbtk's --extension is hardcoded to fna further down.

# Outgroup phylum - must be clearly outside all your genomes' lineages.
# The script looks up the first GTDB representative genome tagged with
# this phylum in your already-downloaded taxonomy file and downloads it
# automatically; you don't need to pick a specific accession by hand.
OUTGROUP_PHYLUM="p__Chlamydiota"

GTDB_RELEASE=release220
GTDB_TAXONOMY_FILE="$PROJECT/gtdbtk_data/${GTDB_RELEASE}/bac120_taxonomy_r220.tsv.gz"
FASTANI_ENV="$PROJECT/envs/fastani_env"   # reusing the datasets CLI installed in step 4
# ---------------------

INDIR=$WORKDIR/gtdbtk_denovo_input
OUTDIR=$WORKDIR/gtdbtk_denovo_out
rm -rf "$INDIR"
mkdir -p "$INDIR" "$OUTDIR"

module load gtdbtk
if [ -f "$PROJECT/gtdbtk_data/env.sh" ]; then
  source "$PROJECT/gtdbtk_data/env.sh"
fi

# --- assemble the input genome set: 6 queries + every genome in $REF_DIR ---
# Extension-agnostic for both: your query genomes are directly in
# $GENOME_DIR (not a genomes/ subfolder) and your reference genomes are
# .fasta not .fna - symlinks are always renamed to "<stem>.fna" in $INDIR
# regardless of the source file's real extension. GTDB-Tk's --extension
# flag only scans ONE extension per run, so mixed extensions on disk
# would otherwise get silently dropped.
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
  # "ref_" prefix is required, not cosmetic: your reference_genomes/ folder
  # is full of real NCBI RefSeq/GenBank accessions (GCF_.../GCA_...) that
  # ARE GTDB's own internal reference genome IDs. Even with --skip_gtdb_refs,
  # gtdbtk still checks input genome IDs against its reference ID list and
  # refuses to run ("same id as GTDB-Tk reference genomes") if they collide.
  ln -sf "$SRC" "$INDIR/ref_${STEM}.fna"
done
echo "Linked ${#REF_FILES[@]} reference genomes from $REF_DIR (id-prefixed 'ref_' to avoid collisions with GTDB-Tk's internal reference genome IDs)."

# --- pick + download one outgroup genome from the GTDB taxonomy file ---
# NOTE: this file is NOT part of the ~100GB package 01_download_gtdbtk_db.slurm
# downloads (that's the internal HMM/pplacer reference data GTDB-Tk itself
# uses) - it's a separate, small (~tens of MB) standalone taxonomy table
# from a different spot on GTDB's server. Download it here if not already
# present rather than assuming step 1 provided it (it doesn't).
if [ ! -f "$GTDB_TAXONOMY_FILE" ]; then
  echo "Taxonomy file not found locally, downloading it (small, not part of the main DB)..."
  mkdir -p "$(dirname "$GTDB_TAXONOMY_FILE")"
  GTDB_RELEASE_SUBDIR=220.0
  TAX_URL="https://data.gtdb.ecogenomic.org/releases/${GTDB_RELEASE}/${GTDB_RELEASE_SUBDIR}/bac120_taxonomy_r220.tsv.gz"
  wget -c "$TAX_URL" -O "$GTDB_TAXONOMY_FILE"
  if ! gzip -t "$GTDB_TAXONOMY_FILE"; then
    echo "ERROR: downloaded taxonomy file failed integrity check (bad URL/404?)." >&2
    echo "Browse https://data.gtdb.ecogenomic.org/releases/${GTDB_RELEASE}/${GTDB_RELEASE_SUBDIR}/ by hand to find the right filename." >&2
    rm -f "$GTDB_TAXONOMY_FILE"
    exit 1
  fi
fi

# NOTE: deliberately NOT using awk's own `exit` after the first match here.
# Under `set -o pipefail`, an early `exit` in awk closes its stdin while
# zcat is still mid-decompression - zcat then gets killed by SIGPIPE (exit
# 141), pipefail propagates that as the pipeline's exit status, and set -e
# kills the whole script instantly with no error message. Letting awk read
# to the end of the (small, tens-of-MB) stream and just skip printing after
# the first match avoids the SIGPIPE entirely; a few extra seconds of zcat
# is a non-issue.
OUTGROUP_LINE=$(zcat "$GTDB_TAXONOMY_FILE" | awk -F'\t' -v p="$OUTGROUP_PHYLUM" '$2 ~ p && !found {print; found=1}')
if [ -z "$OUTGROUP_LINE" ]; then
  echo "ERROR: no genome tagged '$OUTGROUP_PHYLUM' found in $GTDB_TAXONOMY_FILE" >&2
  exit 1
fi
OUTGROUP_ACC_RAW=$(echo "$OUTGROUP_LINE" | cut -f1)
OUTGROUP_ACC=${OUTGROUP_ACC_RAW#RS_}
OUTGROUP_ACC=${OUTGROUP_ACC#GB_}
OUTGROUP_TAXONOMY=$(echo "$OUTGROUP_LINE" | cut -f2)
echo "Outgroup genome: $OUTGROUP_ACC ($OUTGROUP_TAXONOMY)"

# Same ID-collision issue as the reference genomes above applies here too -
# OUTGROUP_ACC was pulled straight out of GTDB's own taxonomy file, so it's
# guaranteed to already be a GTDB reference genome ID. Give the file gtdbtk
# actually scans a prefixed name; OUTGROUP_ACC itself is still used unprefixed
# for the NCBI download/zip/extract paths below (those get cleaned up before
# gtdbtk ever scans $INDIR, so no collision risk there).
OUTGROUP_ID="outgroup_${OUTGROUP_ACC}"
OUTGROUP_FASTA="$INDIR/${OUTGROUP_ID}.fna"
if [ ! -f "$OUTGROUP_FASTA" ]; then
  "$FASTANI_ENV/bin/datasets" download genome accession "$OUTGROUP_ACC" \
    --include genome --filename "$INDIR/${OUTGROUP_ACC}.zip"
  unzip -o -q "$INDIR/${OUTGROUP_ACC}.zip" -d "$INDIR/${OUTGROUP_ACC}_extract"
  # Match whatever extension NCBI's datasets CLI actually used (usually
  # .fna, but don't assume) rather than silently copying nothing.
  # Using find's own -quit (GNU findutils, standard on Ceres) instead of
  # piping into `head -n1` - same SIGPIPE/pipefail risk as the awk fix
  # above if find is still writing when head closes the pipe after one line.
  EXTRACTED=$(find "$INDIR/${OUTGROUP_ACC}_extract" \( -name "*.fna" -o -name "*.fasta" \) -print -quit)
  if [ -z "$EXTRACTED" ]; then
    echo "ERROR: couldn't find a downloaded FASTA under $INDIR/${OUTGROUP_ACC}_extract" >&2
    exit 1
  fi
  cp "$EXTRACTED" "$OUTGROUP_FASTA"
  rm -rf "$INDIR/${OUTGROUP_ACC}.zip" "$INDIR/${OUTGROUP_ACC}_extract"
fi
if [ ! -f "$OUTGROUP_FASTA" ]; then
  echo "ERROR: outgroup genome still missing at $OUTGROUP_FASTA - aborting." >&2
  exit 1
fi

# custom_taxonomy_file: gtdbtk needs this to know which input genome is
# the outgroup when --skip_gtdb_refs means none of these are recognized
# from the GTDB database automatically. Format confirmed from GTDB-Tk's
# own docs: <genome_id (matches filename stem)>\t<full 7-rank taxonomy>.
CUSTOM_TAX_FILE="$INDIR/custom_taxonomy.tsv"
printf '%s\t%s\n' "$OUTGROUP_ID" "$OUTGROUP_TAXONOMY" > "$CUSTOM_TAX_FILE"

# gtdbtk's decorate step (the last stage of de_novo_wf) crashes with
# "IndexError: list index out of range" in _leaf_taxa if ANY leaf has a
# completely empty taxonomy string, instead of just labeling it
# unclassified. Since --skip_gtdb_refs plus the ref_/outgroup_ ID prefixes
# above mean none of your 41 non-outgroup genomes match a known GTDB genome
# ID anymore, every one of them would hit exactly that empty-string case.
# Give each a minimal domain-only placeholder so every leaf has at least
# one non-empty rank - this does NOT assert any real classification, it
# only stops the crash; species/genus for your samples still comes from
# gtdbtk.bac120.summary.tsv (step 3), not from this tree's decoration.
shopt -s nullglob
for F in "$INDIR"/*.fna; do
  GID=$(basename "$F")
  GID=${GID%.fna}
  if [ "$GID" != "$OUTGROUP_ID" ]; then
    printf '%s\td__Bacteria;p__;c__;o__;f__;g__;s__\n' "$GID"
  fi
done >> "$CUSTOM_TAX_FILE"
shopt -u nullglob

echo ""
echo "Genomes going into the de novo tree ($(find "$INDIR" -name '*.fna' | wc -l) total):"
find "$INDIR" -name '*.fna' -exec basename {} \;
echo ""

gtdbtk de_novo_wf \
  --genome_dir "$INDIR" \
  --extension fna \
  --bacteria \
  --outgroup_taxon "$OUTGROUP_PHYLUM" \
  --skip_gtdb_refs \
  --custom_taxonomy_file "$CUSTOM_TAX_FILE" \
  --out_dir "$OUTDIR" \
  --cpus "$SLURM_CPUS_PER_TASK"

echo ""
echo "Tree (Newick): look under $OUTDIR for gtdbtk.bac120.unrooted.tree /"
echo "gtdbtk.bac120.decorated.tree - download and open in iTOL (https://itol.embl.de/)."
echo "Everything except the outgroup will show up unclassified (no species"
echo "label) since we skipped the GTDB reference taxonomy on purpose - the"
echo "labels you care about (your samples + the reference accessions) are"
echo "ones you already know from gtdbtk.bac120.summary.tsv; relabel tips in"
echo "iTOL or with gtdbtk's convert_to_itol if you want nicer display names."
