#!/bin/bash

PROJ_DIR="$HOME/landrace10"
PGGB_DIR="$PROJ_DIR/04_pggb"
MERGE_FA="$PROJ_DIR/02_merged/all_prefixed.fa.gz"
LOG_DIR="$PROJ_DIR/00_logs"

PART_DIR="$PGGB_DIR/rerun_output_all_rice"
mkdir -p "$PART_DIR"

# ------------------------
# Using the default
# ------------------------
PGGB_P="90"
PGGB_S="50000"

echo "Running partition-before-pggb..."
partition-before-pggb -i "$MERGE_FA" -o "$PART_DIR" -n 9 -t 16 -p "$PGGB_P" -s "$PGGB_S" -V 'IRGSP:100000' > "$LOG_DIR/partition_before_pggb_p60.log" 2>&1


# Find all community FASTA files
for part in "$PART_DIR"/*.community.*.fa; do
    # Extract the file basename for the output
    base=$(basename "$part")

    echo "Running PGGB on $base..."
		        
    pggb -i "$part" \
         -o "$PART_DIR/${base}.out" \
         -s "$PGGB_S" -l 250000 -p "$PGGB_P" -c 1 -K 19 -F 0.001 -g 30 \
         -k 23 -f 0 -B 10M \
         -n 9 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ \
         -Y "#" -V IRGSP:100000 --threads 16 --poa-threads 16 \
	 > "$LOG_DIR/pggb_p60.log" 2>&1
done
