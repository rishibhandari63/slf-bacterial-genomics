#!/bin/bash
#SBATCH --job-name=mge_defense
#SBATCH -A YOUR_ACCOUNT_HERE
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o mge_defense_%j.out
#SBATCH -e mge_defense_%j.err

set -euo pipefail

# Mobile genetic elements + defense systems for your 6 isolates:
#   geNomad        - prophages/viruses + plasmid contigs
#   MOB-suite      - plasmid reconstruction/typing (rep/relaxase/mobility)
#   AMRFinderPlus  - AMR + stress/virulence genes (--plus)
#   DefenseFinder  - anti-phage defense systems (RM, Abi, CBASS, etc.)
#   CRISPRCasTyper - CRISPR arrays + cas operon subtyping

# ---- EDIT THESE ----
PROJECT=/project/YOUR_PROJECT_HERE
GENOME_DIR=$PROJECT                 # your 6 sample FASTAs, directly here
PROKKA_DIR=$PROJECT/prokka_out      # for DefenseFinder protein input (preferred)
# --------------------

ENVS=$PROJECT/envs
DBS=$PROJECT/functional_dbs
OUTDIR=$PROJECT/mge_defense_out
mkdir -p "$OUTDIR"

shopt -s nullglob
FA_FILES=("$GENOME_DIR"/*.fasta "$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fa)
shopt -u nullglob
if [ "${#FA_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no genome FASTAs found directly in $GENOME_DIR." >&2
  exit 1
fi
echo "Found ${#FA_FILES[@]} genomes."

# Run each tool so that a failure is RECORDED but doesn't abort the job.
# Rationale from experience on this project: with 5 tools x 6 genomes,
# `set -e` killing the entire run because one tool tripped (AMRFinderPlus
# with a missing -d) threw away the other 4 tools x 6 samples, which had
# nothing wrong with them. Failures are collected and reported in a
# summary at the end, and the job exits nonzero if anything failed - so
# nothing fails silently, it just fails in isolation.
# Note: run_tool always returns 0 so it's safe under `set -e`; the
# FAILURES array is what carries the status.
FAILURES=()
run_tool () {
  local LABEL=$1; shift
  echo "== $LABEL =="
  if ! "$@"; then
    echo "!! FAILED: $LABEL - continuing with remaining tools/samples." >&2
    FAILURES+=("$LABEL")
  fi
  return 0
}

for FA in "${FA_FILES[@]}"; do
  NAME=$(basename "$FA"); NAME=${NAME%.*}
  S=$OUTDIR/$NAME
  mkdir -p "$S"

  # --- geNomad (prophages + plasmid classification) ---
  if [ -d "$S/genomad" ] && [ -n "$(find "$S/genomad" -name '*_summary' -o -name '*summary*' 2>/dev/null | head -c1)" ]; then
    echo "== geNomad: $NAME already done, skipping =="
  else
    rm -rf "$S/genomad"
    run_tool "geNomad: $NAME" \
      env PATH="$ENVS/mge_env/bin:$PATH" "$ENVS/mge_env/bin/genomad" end-to-end --cleanup \
        --threads "$SLURM_CPUS_PER_TASK" \
        "$FA" "$S/genomad" "$DBS/genomad/genomad_db"
  fi

  # --- MOB-suite (plasmid reconstruction/typing) ---
  if [ -f "$S/mob_recon/contig_report.txt" ]; then
    echo "== MOB-suite: $NAME already done, skipping =="
  else
    rm -rf "$S/mob_recon"
    run_tool "MOB-suite: $NAME" \
      env PATH="$ENVS/mge_env/bin:$PATH" "$ENVS/mge_env/bin/mob_recon" \
        -i "$FA" -o "$S/mob_recon" -n "$SLURM_CPUS_PER_TASK" --force
  fi

  # --- AMRFinderPlus (AMR/stress/virulence) ---
  if [ -s "$S/amrfinder.tsv" ]; then
    echo "== AMRFinderPlus: $NAME already done, skipping =="
  else
    # -d points at the explicitly-located database from script 12 ("latest"
    # symlink is maintained by amrfinder_update). Without it, amrfinder
    # falls back to a nonexistent build-time path when $CONDA_PREFIX is
    # unset (same failure we fixed in the download job). If you see
    # "Reverting to hard coded directory: /opt/conda/..." in the log, the
    # -d below did not reach the binary - you're running an older copy of
    # this script.
    run_tool "AMRFinderPlus: $NAME" \
      env PATH="$ENVS/mge_env/bin:$PATH" "$ENVS/mge_env/bin/amrfinder" \
        -n "$FA" --plus --threads "$SLURM_CPUS_PER_TASK" \
        -d "$DBS/amrfinderplus/latest" \
        -o "$S/amrfinder.tsv"
  fi

  # --- DefenseFinder (anti-phage systems) ---
  # Prefer Prokka's .faa (proteins in genomic order, which DefenseFinder
  # needs); fall back to the nucleotide genome, which it also accepts.
  DF_IN=""
  shopt -s nullglob
  CAND=("$PROKKA_DIR/$NAME"/*.faa "$PROKKA_DIR/${NAME}.faa")
  shopt -u nullglob
  [ "${#CAND[@]}" -gt 0 ] && DF_IN=${CAND[0]} || DF_IN=$FA
  if [ -d "$S/defensefinder" ] && [ -n "$(ls -A "$S/defensefinder" 2>/dev/null)" ]; then
    echo "== DefenseFinder: $NAME already done, skipping =="
  else
    rm -rf "$S/defensefinder"
    mkdir -p "$S/defensefinder"
    run_tool "DefenseFinder: $NAME (input: $DF_IN)" \
      env PATH="$ENVS/defense_env/bin:$PATH" "$ENVS/defense_env/bin/defense-finder" run \
        -o "$S/defensefinder" \
        --models-dir "$DBS/defensefinder_models" \
        -w "$SLURM_CPUS_PER_TASK" \
        "$DF_IN"
  fi

  # --- CRISPRCasTyper ---
  if [ -d "$S/cctyper" ] && [ -n "$(ls -A "$S/cctyper" 2>/dev/null)" ]; then
    echo "== cctyper: $NAME already done, skipping =="
  else
    rm -rf "$S/cctyper"   # cctyper requires its output dir not to exist
    run_tool "cctyper: $NAME" \
      env PATH="$ENVS/defense_env/bin:$PATH" "$ENVS/defense_env/bin/cctyper" \
        "$FA" "$S/cctyper" -t "$SLURM_CPUS_PER_TASK"
  fi
done

echo ""
echo "Done. Per sample under $OUTDIR/<sample>/:"
echo "  genomad/         - *_virus_summary.tsv (prophages), *_plasmid_summary.tsv"
echo "  mob_recon/       - contig_report.txt (chromosome vs plasmid), mobtyper results"
echo "  amrfinder.tsv    - AMR/stress/virulence hits (empty-ish is a fine result for"
echo "                     environmental isolates - report it either way)"
echo "  defensefinder/   - defense_finder_systems.tsv (one row per detected system)"
echo "  cctyper/         - CRISPR_Cas.tab + crisprs_all.tab"
echo ""
echo "Note: DefenseFinder's -w flag = workers; if your installed version"
echo "rejects it, drop that flag (it then uses its default parallelism)."

echo ""
if [ "${#FAILURES[@]}" -eq 0 ]; then
  echo "All tools completed for all genomes."
else
  echo "!! ${#FAILURES[@]} tool run(s) FAILED (everything else completed):" >&2
  printf '   - %s\n' "${FAILURES[@]}" >&2
  echo "" >&2
  echo "Fix the cause and resubmit - completed tools/samples are skip-detected," >&2
  echo "so only the failed ones will rerun." >&2
  exit 1
fi
