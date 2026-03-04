#!/usr/bin/env bash

# Usage
usage() {
    echo "Usage: run_panacus.sh -s SAMPLE_ORDER [-t THREADS]"
    echo
    echo "  -s SAMPLE_ORDER   Sample order file (required)"
    echo "  -t THREADS        Number of threads (optional, default=2)"
    exit 1
}

THREADS=2
SAMPLE_ORDER=""

# Parse command line arguments
while getopts "s:t:" opt; do
    case $opt in
        s) SAMPLE_ORDER="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check if supplied
if [ -z "$SAMPLE_ORDER" ]; then
    echo "Error: SAMPLE_ORDER is required."
    usage
fi

# Check sample order file
if [ ! -f "$SAMPLE_ORDER" ]; then    
    echo "Error: Sample order file '$SAMPLE_ORDER' does not exist."
    exit 1
fi

# Process all output directories
for DIR in output_all_Omer_rerun_p*; do
    echo "Processing $DIR..."
    cd "$DIR"

    # 1. find the GFA file
    GFA=$(ls *.gfa)
    PREFIX=$(basename "$GFA" .gfa)

    # 2. graph statistics using panacus
    panacus info "$GFA" -S -t "$THREADS" > "${PREFIX}.info.tsv"

    # 3a. ordered histgrowth (bp)
    RUST_LOG=info panacus ordered-histgrowth \
        -c bp \
        -t "$THREADS" \
        -l 1,2,7 \
        -S \
        -O "../$SAMPLE_ORDER" \
        "$GFA" \
        > "${PREFIX}.ordered-histgrowth.bp.tsv"

    # 3b. ordered histgrowth (node)
    RUST_LOG=info panacus ordered-histgrowth \
        -c node \
        -t "$THREADS" \
        -l 1,2,7 \
        -S \
        -O "../$SAMPLE_ORDER" \
        "$GFA" \
        > "${PREFIX}.ordered-histgrowth.node.tsv"
		
	# 4. visualize
    panacus-visualize \
	"${PREFIX}.ordered-histgrowth.bp.tsv" \
	> "${PREFIX}.ordered-histgrowth.bp.pdf"
	
	# 4. visualize
    panacus-visualize \
	"${PREFIX}.ordered-histgrowth.node.tsv" \
	> "${PREFIX}.ordered-histgrowth.node.pdf"
    
    # Exit the directory
    cd ..
done

echo "All Panacus-related steps succesful..."
