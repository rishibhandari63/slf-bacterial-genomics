#!/bin/bash
#SBATCH --job-name=funct_dbs
#SBATCH -A YOUR_ACCOUNT_HERE
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o funct_dbs_%j.out
#SBATCH -e funct_dbs_%j.err

set -euo pipefail

# One-time download of all databases for the functional-genomics extension.
# Sizes (approximate, checked Aug 2026): eggNOG ~50 GB, dbCAN ~10 GB,
# antiSMASH ~9 GB, geNomad ~1.5 GB, AMRFinderPlus ~1 GB, DefenseFinder
# models ~100 MB, MOB-suite ~3 GB. Total ~75 GB - make sure your /project
# quota has room (my_quota) BEFORE submitting.
#
# Every step below is idempotent - safe to resubmit if the job times out
# partway (already-complete downloads are detected and skipped).

# ---- EDIT THIS ----
PROJECT=/project/YOUR_PROJECT_HERE
# --------------------

ENVS=$PROJECT/envs
DBS=$PROJECT/functional_dbs
mkdir -p "$DBS"

# --- eggNOG-mapper data (~50 GB) ---
EGGNOG_DIR=$DBS/eggnog
mkdir -p "$EGGNOG_DIR"
if [ -f "$EGGNOG_DIR/eggnog.db" ] && [ -f "$EGGNOG_DIR/eggnog_proteins.dmnd" ]; then
  echo "== eggNOG data already present, skipping =="
else
  echo "== downloading eggNOG data (largest download - ~50 GB) =="
  # The download host hardcoded in eggnog-mapper's downloader moved from
  # eggnogdb.embl.de (dead - DNS no longer resolves) to eggnog5.embl.de;
  # patch the installed script in place before running it. Idempotent:
  # sed is a no-op once patched. (Known upstream issue, eggnog-mapper #589.)
  sed -i 's|eggnogdb\.embl\.de|eggnog5.embl.de|g' \
    "$ENVS/funct_env/bin/download_eggnog_data.py"
  "$ENVS/funct_env/bin/download_eggnog_data.py" -y --data_dir "$EGGNOG_DIR"
fi

# --- dbCAN v5 databases (~10 GB) ---
DBCAN_DIR=$DBS/dbcan
if [ -d "$DBCAN_DIR" ] && [ -n "$(ls -A "$DBCAN_DIR" 2>/dev/null)" ]; then
  echo "== dbCAN db dir non-empty, assuming present (delete $DBCAN_DIR to force redownload) =="
else
  echo "== downloading dbCAN databases =="
  # --aws_s3 is the documented faster/more stable mirror; --resume makes
  # interrupted downloads restartable.
  "$ENVS/funct_env/bin/run_dbcan" database --db_dir "$DBCAN_DIR" --aws_s3 --resume
fi

# --- antiSMASH databases (~9 GB) ---
ANTISMASH_DIR=$DBS/antismash
mkdir -p "$ANTISMASH_DIR"
# NOTE: deliberately NO "skip if directory non-empty" check here. An earlier
# version of this script had one, and it caused a real failure: the first
# run died partway (hmmpress not on PATH) leaving a PARTIALLY populated
# directory, then every rerun saw "non-empty" and skipped the downloader -
# so Pfam was never pressed and antiSMASH 8's MITE dataset was never
# fetched at all, which only surfaced later as
# "ValueError: No matching database in location .../antismash/mite".
# download-antismash-databases does its own per-file presence/checksum
# checks and skips what's already complete, so just always run it.
echo "== downloading/verifying antiSMASH databases =="
# It shells out to hmmpress (same env's bin), so the env bin must be on
# PATH - same subprocess-PATH lesson as pyani/nucmer.
PATH="$ENVS/antismash_env/bin:$PATH" \
  "$ENVS/antismash_env/bin/download-antismash-databases" --database-dir "$ANTISMASH_DIR"

# Post-check: antiSMASH 8 refuses to start if any expected dataset is
# missing, and the failure appears at RUN time (script 14), not here -
# so verify the ones that have bitten us before, while we can still act.
for SUB in pfam mite; do
  if [ ! -d "$ANTISMASH_DIR/$SUB" ] || [ -z "$(ls -A "$ANTISMASH_DIR/$SUB" 2>/dev/null)" ]; then
    echo "WARNING: $ANTISMASH_DIR/$SUB is missing or empty after the download." >&2
    echo "  If script 14 later fails with 'No matching database in location'," >&2
    echo "  check whether your antiSMASH version ($("$ENVS/antismash_env/bin/antismash" --version 2>&1)) needs" >&2
    echo "  a newer downloader, and see https://docs.antismash.secondarymetabolites.org/install/" >&2
  fi
done

# --- geNomad database (~1.5 GB) ---
GENOMAD_DIR=$DBS/genomad
if [ -d "$GENOMAD_DIR/genomad_db" ]; then
  echo "== geNomad db already present, skipping =="
else
  echo "== downloading geNomad database =="
  mkdir -p "$GENOMAD_DIR"
  "$ENVS/mge_env/bin/genomad" download-database "$GENOMAD_DIR"
fi

# --- AMRFinderPlus database (~1 GB; versioned, updates in place) ---
# Called by full path with no conda activation, amrfinder finds $CONDA_PREFIX
# unset and falls back to a directory hardcoded at bioconda BUILD time
# (/opt/conda/conda-bld/...) which doesn't exist here - so give it an
# explicit database directory instead of using `amrfinder -u`. Script 15
# passes the same directory at runtime via -d.
AMRF_DB=$DBS/amrfinderplus
mkdir -p "$AMRF_DB"
echo "== updating AMRFinderPlus database =="
"$ENVS/mge_env/bin/amrfinder_update" -d "$AMRF_DB"
# Note (informational): NCBI's newest AMR database (2026-08+) requires
# amrfinder >=4.2; the conda-installed 3.12.x fetches the newest COMPATIBLE
# database version instead, which is fine to use and cite. To get the very
# latest DB, upgrade the tool first:
#   mamba install -p $ENVS/mge_env 'ncbi-amrfinderplus>=4.2'   (via salloc)

# --- DefenseFinder models (~100 MB) ---
DF_MODELS=$DBS/defensefinder_models
mkdir -p "$DF_MODELS"
echo "== updating DefenseFinder models =="
# macsypy (under defense-finder) checks for hmmsearch on PATH even during
# a model update - hmmsearch lives in the same env's bin, so prefix PATH
# (same lesson as antiSMASH/hmmpress above).
PATH="$ENVS/defense_env/bin:$PATH" \
  "$ENVS/defense_env/bin/defense-finder" update --models-dir "$DF_MODELS"

# --- MOB-suite databases (~3 GB; mob_init fetches to the env's data dir) ---
echo "== initializing MOB-suite databases =="
PATH="$ENVS/mge_env/bin:$PATH" "$ENVS/mge_env/bin/mob_init" || echo "mob_init returned nonzero - if mob_recon later complains about databases, rerun this step."

echo ""
echo "All databases done. Layout:"
du -sh "$DBS"/* 2>/dev/null || true
echo ""
echo "Next: sbatch 13_funct_annotation.slurm (after Prokka .faa files exist),"
echo "      sbatch 14_antismash.slurm, sbatch 15_mge_defense.slurm"
