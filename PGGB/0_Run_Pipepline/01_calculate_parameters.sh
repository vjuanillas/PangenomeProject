#!/bin/bash

# Directories
PROJ_DIR="$HOME/landrace10"
MERGE_FA="$PROJ_DIR/02_merged/all_prefixed.fa.gz"
REPEAT_FILE="/opt/home/venice/landrace10_asm5/repeatmasker/all_rice_O_mer.fa.out" # specify RepeatMasker out file, not tbl
OUT_DIR="$PROJ_DIR/03_mash"

echo "Running Mash triangle (default)..."
mash triangle -i "$MERGE_FA" > "$OUT_DIR/all_rice.mash_triangle.default.txt"

echo "Running Mash triangle with -s 10000..."
mash triangle -s 10000 -i "$MERGE_FA" > "$OUT_DIR/all_rice.mash_triangle.10k.txt"

echo "Top two divergence values (default sketch):"
sed 1,1d "$OUT_DIR/all_rice.mash_triangle.default.txt" | tr '\t' '\n' | grep -v chr | LC_ALL=C sort -g -k1nr | uniq | head -n 2

echo "Top two divergence values (sketch=10000):"
sed 1,1d "$OUT_DIR/all_rice.mash_triangle.10k.txt" | tr '\t' '\n' | grep -v chr | LC_ALL=C sort -g -k1nr | uniq | head -n 2

# ----------------------------
# Suggested PGGB -p based on maximum divergence
# ----------------------------
max_dist=$(sed 1,1d "$OUT_DIR/all_rice.mash_triangle.10k.txt" | tr '\t' '\n' | grep -v chr | LC_ALL=C sort -g -k1nr | uniq | head -n 1)

# Convert Mash distance to % identity
PGGB_P=$(awk -v d="$max_dist" 'BEGIN {printf "%d", (1-d)*100}')

echo "Suggested PGGB -p: $PGGB_P"

# ----------------------------
# Suggested PGGB -s based on maximum repeat length
# ----------------------------

# Compute max repeat length based on Repeatmasker output file
max_repeat=$(sed 1,3d "$REPEAT_FILE" | awk '{len=$7-$6; if(len>max) max=len} END{print max}')
echo "Max repeat length: $max_repeat"

# Add 1kb buffer for PGGB
PGGB_S=$((max_repeat + 1000))
echo "Suggested PGGB -s: $PGGB_S"
