#!/bin/bash
# One-time, INTERACTIVE setup for pyani (v0.2.x - the stable/archived release;
# the exact tool that produces the clustered-heatmap-with-dendrogram figure
# style you're after). Conda must run on a compute node, not the login node:
#
#   1. From the Ceres login node:
#        salloc -A YOUR_ACCOUNT_HERE -p short -c 4 --mem=16G -t 01:00:00
#   2. Once on a compute node:
#        bash 08_install_pyani_env.sh
#
# IMPORTANT: pyani's ANIm method is pinned to MUMmer 3.23 specifically -
# the pyani docs say it has "not yet been tested with MUMmer 4.x". Ceres's
# `module load mummer` gives you 4.0.0rc1, which is NOT what we want here,
# so this installs its own MUMmer 3.23 inside the env instead of touching
# that module at all (also sidesteps mixing Lmod modules with conda PATH
# handling, which has already bitten us once in this project).

set -euo pipefail

# ---- EDIT THIS ----
PROJECT=/project/YOUR_PROJECT_HERE
# --------------------

ENV_PREFIX=$PROJECT/envs/pyani_env

module load miniconda

if [ -d "$ENV_PREFIX" ]; then
  echo "Environment already exists at $ENV_PREFIX, skipping creation."
else
  mkdir -p "$PROJECT/envs"
  # python=3.8 because pyani 0.2.x's own docs pin to this generation of
  # Python; mummer=3.23 per pyani's documented ANIm requirement above.
  mamba create -y --prefix "$ENV_PREFIX" -c bioconda -c conda-forge \
    python=3.8 pyani "mummer=3.23"
fi

echo ""
echo "Checking the install (calling binaries directly by full path - source"
echo "activate has been unreliable on this cluster in this session, so"
echo "09_run_pyani_anim.slurm does the same):"
"$ENV_PREFIX/bin/average_nucleotide_identity.py" --version || true
"$ENV_PREFIX/bin/nucmer" --version || true

echo ""
echo "Done. 09_run_pyani_anim.slurm will call $ENV_PREFIX/bin/average_nucleotide_identity.py directly."
