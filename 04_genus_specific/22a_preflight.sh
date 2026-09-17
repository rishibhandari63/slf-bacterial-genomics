#!/bin/bash
# Read-only survey of the Massilia directory. Run it on the login node -
# it takes seconds and writes nothing.
#
#   bash 22a_preflight.sh
#
# It tells you exactly what to put in REF_DIR, and whether Sample4's own
# proteome is sitting in the reference set (which would silently produce a
# zero result).

PROJECT=/project/kanglab/rishi.bhandari
GENUS_DIR=$PROJECT/massilia
SAMPLE=Sample4
ISOLATE_FAA=$PROJECT/PROKKA_DIR/prokka_out/$SAMPLE/$SAMPLE.faa

echo "=================================================================="
echo " Pre-flight for script 22  ($SAMPLE vs $GENUS_DIR)"
echo "=================================================================="

# ---- 1. the isolate proteome the eggNOG join depends on -----------------
echo
echo "[1] Isolate proteome (must be the one eggNOG-mapper annotated)"
if [ -s "$ISOLATE_FAA" ]; then
    FIRST_ID=$(grep -m1 '^>' "$ISOLATE_FAA" | cut -d' ' -f1)
    echo "    OK   $ISOLATE_FAA"
    echo "         $(grep -c '^>' "$ISOLATE_FAA") proteins, first ID ${FIRST_ID#>}"
else
    echo "    MISSING: $ISOLATE_FAA"
    echo "    Find it with:  find $PROJECT -name '$SAMPLE.faa' 2>/dev/null"
    FIRST_ID=""
fi

EGG=$PROJECT/funct_annotation_out/eggnog/$SAMPLE/$SAMPLE.emapper.annotations
if [ -s "$EGG" ]; then
    EGG_ID=$(awk -F'\t' '!/^#/{print $1; exit}' "$EGG")
    echo "    eggNOG first ID: $EGG_ID"
    if [ -n "$FIRST_ID" ] && [ "${FIRST_ID#>}" != "$EGG_ID" ]; then
        echo "    !! LOCUS TAGS DIFFER between the .faa and the eggNOG table."
        echo "       The functional join will return 0%. These must be the same"
        echo "       Prokka run."
    fi
fi

# ---- 2. where the reference proteomes actually are ----------------------
echo
echo "[2] .faa files under $GENUS_DIR"
if [ ! -d "$GENUS_DIR" ]; then
    echo "    Directory not found."
    exit 1
fi
mapfile -t FAAS < <(find "$GENUS_DIR" -name '*.faa' -type f 2>/dev/null | sort)
echo "    found ${#FAAS[@]} .faa files"

if [ "${#FAAS[@]}" -eq 0 ]; then
    echo
    echo "    No protein FASTA here. Script 22 needs .faa, not .gff."
    echo "    Look for the full Prokka output:"
    echo "      find $GENUS_DIR -name '*.faa' -o -name 'PROKKA*' -type d | head"
    echo "    If you kept only the GFFs, the proteins are still inside them -"
    echo "    ask and I will give you a GFF-to-faa extractor."
    exit 1
fi

echo
echo "    Directories holding them (counts):"
printf '%s\n' "${FAAS[@]}" | xargs -n1 dirname | sort | uniq -c | sort -rn | head -8

# The common parent is what REF_DIR should be: the glob in script 22 is
# recursive, so the top of the tree is the right answer.
COMMON=$(printf '%s\n' "${FAAS[@]}" | xargs -n1 dirname | sed "s|$GENUS_DIR||" \
         | awk -F/ '{print $2}' | sort -u)
NLVL=$(echo "$COMMON" | grep -c .)
if [ "$NLVL" -eq 1 ] && [ -n "$COMMON" ]; then
    SUGGEST="$GENUS_DIR/$COMMON"
else
    SUGGEST="$GENUS_DIR"
fi

# ---- 3. is the isolate inside the reference set? ------------------------
echo
echo "[3] Is $SAMPLE inside its own reference set?"
HITS=0
for F in "${FAAS[@]}"; do
    WHY=""
    case "$(basename "$F")" in "$SAMPLE".faa) WHY="basename" ;; esac
    case "$F" in *"$SAMPLE"*) WHY="${WHY:-path}" ;; esac
    if [ -n "$FIRST_ID" ] && \
       [ "$(grep -m1 '^>' "$F" | cut -d' ' -f1)" = "$FIRST_ID" ]; then
        WHY="${WHY:+$WHY, }locus tags"
    fi
    if [ -n "$WHY" ]; then
        echo "    YES ($WHY): $F"
        HITS=$((HITS + 1))
    fi
done
if [ "$HITS" -eq 0 ]; then
    echo "    No copy of $SAMPLE found in the reference set."
    echo "    Note: it SHOULD be in the Roary GFF set but must NOT be here."
else
    echo "    -> script 22 will exclude these automatically. Confirm the count"
    echo "       above matches what you expect."
fi

# ---- 4. what to do ------------------------------------------------------
USABLE=$(( ${#FAAS[@]} - HITS ))
echo
echo "=================================================================="
echo " Set this in script 22:"
echo "     Sample4) REF_DIR=$SUGGEST ;;"
echo
echo " Comparator genomes after exclusions: $USABLE"
echo " (previously 30 Telluria - cite the new n everywhere in the paper)"
echo "=================================================================="
