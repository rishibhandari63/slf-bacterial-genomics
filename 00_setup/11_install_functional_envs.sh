#!/bin/bash
# One-time, INTERACTIVE setup of conda environments for the functional /
# comparative genomics extension of this pipeline:
#
#   funct_env    - eggNOG-mapper (KEGG/COG functional annotation)
#                  + run_dbcan v5 (CAZymes)
#   antismash_env - antiSMASH (secondary metabolite / BGC detection)
#   mge_env      - geNomad (prophages/plasmids), MOB-suite (plasmid typing),
#                  AMRFinderPlus (AMR/virulence genes)
#   defense_env  - DefenseFinder (anti-phage defense systems)
#                  + CRISPRCasTyper (cctyper; CRISPR-Cas arrays)
#
# Separate envs on purpose: antiSMASH in particular has a heavy, fussy
# dependency stack that regularly conflicts with other bioinformatics
# packages - do not merge these to save space.
#
# As established earlier in this project: conda must run on a COMPUTE node.
#   1. From the Ceres login node:
#        salloc -A YOUR_ACCOUNT_HERE -p short -c 4 --mem=16G -t 02:00:00
#   2. Once on a compute node:
#        bash 11_install_functional_envs.sh
#
# Databases are NOT downloaded here - that's 12_download_functional_dbs.slurm
# (they total tens of GB; a batch job with a real time limit is the right
# place for that, not an interactive session).

set -euo pipefail

# ---- EDIT THIS ----
PROJECT=/project/YOUR_PROJECT_HERE
# --------------------

ENVS=$PROJECT/envs
mkdir -p "$ENVS"

module load miniconda

create_env () {
  local PREFIX=$1; shift
  if [ -d "$PREFIX" ]; then
    echo "== $PREFIX already exists, skipping creation =="
  else
    echo "== creating $PREFIX =="
    mamba create -y --prefix "$PREFIX" -c bioconda -c conda-forge "$@"
  fi
}

# eggNOG-mapper + dbCAN share an env fine (both python, no native conflicts).
create_env "$ENVS/funct_env"     eggnog-mapper dbcan

# antiSMASH pulls its own pinned deps - keep isolated.
create_env "$ENVS/antismash_env" antismash

# MGE toolkit.
create_env "$ENVS/mge_env"       genomad mob_suite ncbi-amrfinderplus

# Defense systems. defense-finder is on bioconda; cctyper ships its own
# CRISPR-Cas typing models.
create_env "$ENVS/defense_env"   defense-finder cctyper

echo ""
echo "== sanity checks (full-path invocation, as everywhere in this pipeline) =="
"$ENVS/funct_env/bin/emapper.py" --version || true
"$ENVS/funct_env/bin/run_dbcan" --help >/dev/null 2>&1 && echo "run_dbcan: OK" || echo "run_dbcan: CHECK INSTALL"
"$ENVS/antismash_env/bin/antismash" --version || true
"$ENVS/mge_env/bin/genomad" --version || true
"$ENVS/mge_env/bin/mob_recon" --version || true
"$ENVS/mge_env/bin/amrfinder" --version || true
"$ENVS/defense_env/bin/defense-finder" --version || true
"$ENVS/defense_env/bin/cctyper" --version || true

echo ""
echo "Done. Next: sbatch 12_download_functional_dbs.slurm"
