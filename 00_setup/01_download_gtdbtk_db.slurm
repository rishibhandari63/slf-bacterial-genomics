#!/bin/bash
#SBATCH --job-name=gtdbtk_db_download
#SBATCH -A YOUR_ACCOUNT_HERE
#SBATCH --partition=medium
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -o gtdbtk_db_download_%j.out
#SBATCH -e gtdbtk_db_download_%j.err

set -euo pipefail

# ---- EDIT THESE ----
PROJECT=/project/YOUR_PROJECT_HERE
# Release must match your loaded gtdbtk module version. As of this writing,
# Ceres's newest gtdbtk module is 2.4.0, which pairs with GTDB release R220
# (compatible range: gtdbtk 2.4.0-2.6.1 <-> R220). Check the compatibility
# table at https://ecogenomics.github.io/GTDBTk/installing/index.html before
# you run this if you're using a different module version.
GTDB_RELEASE=release220
GTDB_RELEASE_SUBDIR=220.0
GTDB_RELEASE_NUM=220
# ---------------------

DB_DIR="$PROJECT/gtdbtk_data/${GTDB_RELEASE}"
mkdir -p "$DB_DIR"
cd "$DB_DIR"

module load gtdbtk
echo "Loaded gtdbtk version:"
gtdbtk -v || true

# Note: for a numbered release (as opposed to "latest"), the archive's
# filename embeds the release number, e.g. gtdbtk_r220_data.tar.gz - it is
# NOT just gtdbtk_data.tar.gz (that name is only used under .../latest/).
URL="https://data.gtdb.ecogenomic.org/releases/${GTDB_RELEASE}/${GTDB_RELEASE_SUBDIR}/auxillary_files/gtdbtk_package/full_package/gtdbtk_r${GTDB_RELEASE_NUM}_data.tar.gz"
echo "Downloading: $URL"

# -c makes this resumable if the job times out partway through and you
# need to resubmit; this file is roughly 100 GB, so budget real time for it.
wget -c "$URL" -O gtdbtk_data.tar.gz

# Sanity check before spending time extracting a partial/corrupt file.
if ! gzip -t gtdbtk_data.tar.gz; then
  echo "ERROR: downloaded archive failed integrity check. If this was a" >&2
  echo "URL 404 or redirect, browse https://data.gtdb.ecogenomic.org/releases/${GTDB_RELEASE}/ " >&2
  echo "by hand (e.g. from a login-node browser tab) to find the correct path" >&2
  echo "and fix GTDB_RELEASE / GTDB_RELEASE_SUBDIR / URL above." >&2
  exit 1
fi

echo "Extracting..."
tar -xzf gtdbtk_data.tar.gz --strip-components=1

echo "Done. Set this in your environment / other job scripts:"
echo "export GTDBTK_DATA_PATH=$DB_DIR"

# Leave a small marker file the other scripts can source.
echo "export GTDBTK_DATA_PATH=$DB_DIR" > "$PROJECT/gtdbtk_data/env.sh"
