#!/bin/bash
#SBATCH --job-name=massilia_roary
#SBATCH -A kanglab
#SBATCH -p ceres
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=2-00:00:00
#SBATCH -o logs/roary_%j.out
#SBATCH -e logs/roary_%j.err

set -euo pipefail

########## config ##########
WORKDIR="/project/kanglab/rishi.bhandari/massilia/core_phylo"
GFFDIR="$WORKDIR/gff"
OUTDIR="$WORKDIR/roary_i80"       # must NOT exist; Roary refuses to overwrite
CORE_IDENT=80                     # BLASTp %id for clustering
CORE_DEF=95                       # % of genomes a gene must be in to count as core
GROUP_LIMIT=600000                # default 50000 is far too low for a genus
ISOLATE="Sample4"                 # must be present in $GFFDIR or the tree is pointless
############################

THREADS="${SLURM_CPUS_PER_TASK:-8}"

conda deactivate 2>/dev/null || true
module purge
module load roary/3.13.0--pl526h516909a_0

# ---- fail fast on a missing MAFFT: with -e --mafft, Roary needs it, and
# ---- discovering that two days into the job is expensive.
command -v mafft >/dev/null 2>&1 || module load mafft 2>/dev/null || true
if ! command -v mafft >/dev/null 2>&1; then
    echo "ERROR: mafft not on PATH. 'roary -e --mafft' will fail." >&2
    exit 1
fi
echo "mafft: $(command -v mafft)"

mkdir -p logs
cd "$WORKDIR"

# ---- the isolate must actually be in the input set ----------------------
if ! ls "$GFFDIR"/*"$ISOLATE"*.gff >/dev/null 2>&1; then
    echo "ERROR: no GFF matching '$ISOLATE' in $GFFDIR." >&2
    echo "       Building a genus tree without your own isolate in it is the" >&2
    echo "       one mistake that wastes the whole run." >&2
    exit 1
fi
echo "isolate GFF: $(ls "$GFFDIR"/*"$ISOLATE"*.gff)"

# ---- every GFF must be Prokka-style: annotation followed by ##FASTA -----
# NCBI .gff files have no sequence block and Roary will produce nonsense.
BAD=0
for G in "$GFFDIR"/*.gff; do
    grep -q '^##FASTA' "$G" || { echo "  NO ##FASTA BLOCK: $G" >&2; BAD=$((BAD+1)); }
done
if [ "$BAD" -gt 0 ]; then
    echo "ERROR: $BAD GFF file(s) carry no sequence block." >&2
    echo "       Roary needs Prokka-style GFF3 (annotation + ##FASTA + contigs)." >&2
    echo "       Re-annotate those genomes with Prokka rather than using the" >&2
    echo "       NCBI .gff directly." >&2
    exit 1
fi

n=$(ls "$GFFDIR"/*.gff | wc -l)
echo ">> Roary on $n genomes | -i $CORE_IDENT -cd $CORE_DEF -g $GROUP_LIMIT"

roary -e --mafft \
    -p "$THREADS" \
    -i "$CORE_IDENT" \
    -cd "$CORE_DEF" \
    -g "$GROUP_LIMIT" \
    -f "$OUTDIR" \
    "$GFFDIR"/*.gff

echo "--- summary ---"
cat "$OUTDIR/summary_statistics.txt" || true

ALN="$OUTDIR/core_gene_alignment.aln"
if [ ! -s "$ALN" ]; then
    echo ">> No alignment produced. Check $OUTDIR/summary_statistics.txt"
    exit 1
fi

taxa=$(grep -c '^>' "$ALN")
cols=$(awk '/^>/{if(s){print length(s); exit}} !/^>/{s=s $0} END{if(s)print length(s)}' "$ALN")
core=$(awk -F'\t' '/^Core genes/{print $NF}' "$OUTDIR/summary_statistics.txt" || echo "?")

echo ">> SUCCESS: $ALN"
ls -lh "$ALN"
echo "   Taxa in alignment : $taxa  (expected $n)"
echo "   Alignment columns : $cols"
echo "   Core genes        : $core"

# NOTE: plain '[ x -ne y ] && echo ...' as the final statement makes the whole
# script exit 1 whenever the counts DO match, and Slurm then reports the job as
# FAILED after a successful two-day run. Use a full if-block.
if [ "$taxa" -ne "$n" ]; then
    echo "   WARNING: taxon count does not match genome count"
fi

# ---- did the clustering actually work this time? -----------------------
if [ "$core" != "?" ] && [ "$core" -lt 200 ] 2>/dev/null; then
    echo ""
    echo "   !! Only $core core genes from $n genomes. At -i 95 this dataset gave"
    echo "      89. If -i 80 has not moved it well into the hundreds, Roary is"
    echo "      still over-splitting and the alignment is not a sound basis for"
    echo "      a tree. Switch to a marker-gene approach (GToTree, or GTDB-Tk"
    echo "      de_novo_wf on the bac120 set) rather than lowering -i further."
fi

echo ""
echo "Next: RAxML on $ALN, then compare the Massilia clade order against the"
echo "bac120 tree in Figure 1b. Two independent methods agreeing is the"
echo "cheapest validation available for this tree."
exit 0
