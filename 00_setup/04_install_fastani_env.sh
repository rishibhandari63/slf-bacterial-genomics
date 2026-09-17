#!/bin/bash
# One-time, INTERACTIVE setup for a small conda environment with fastANI,
# the NCBI datasets CLI (for pulling reference genomes to compare against),
# and pandas/matplotlib (for turning fastANI's output into a heatmap via
# 07_plot_ani_heatmap.py). Conda must run on a compute node, not the login
# node, so:
#
#   1. From the Ceres login node, request an interactive session, e.g.:
#        salloc -A YOUR_ACCOUNT_HERE -p short -c 4 --mem=16G -t 01:00:00
#   2. Once it drops you onto a compute node, run this script:
#        bash 04_install_fastani_env.sh
#
# This only needs to be done once; the resulting env is reused by
# 05_download_refs_and_fastani.slurm.

set -euo pipefail

# ---- EDIT THIS ----
PROJECT=/project/YOUR_PROJECT_HERE
# --------------------

ENV_PREFIX=$PROJECT/envs/fastani_env

module load miniconda

if [ -d "$ENV_PREFIX" ]; then
  echo "Environment already exists at $ENV_PREFIX."
  echo "Making sure pandas/matplotlib are in it too (no-op if already present)..."
  mamba install -y --prefix "$ENV_PREFIX" -c bioconda -c conda-forge \
    pandas matplotlib
else
  mkdir -p "$PROJECT/envs"
  mamba create -y --prefix "$ENV_PREFIX" -c bioconda -c conda-forge \
    fastani ncbi-datasets-cli pandas matplotlib
fi

echo ""
echo "Done. Activate it in job scripts with:"
echo "  module load miniconda"
echo "  source activate $ENV_PREFIX"
