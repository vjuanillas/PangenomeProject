#!/bin/bash

# activate LTR_retriever
# replce the paths to tools if necessary

# Ensure the script exits on any error
set -e

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <assembly_file> <index_name> <threads>"
    exit 1
fi

# Read input arguments
ASSEMBLY_PATH=$1
INDEX_NAME=$2
THREADS=$3

# Extract the base name of the assembly file (without path and extension)
ASSEMBLY_BASENAME=$(basename "$ASSEMBLY_PATH" .fasta)

# Define paths to tools
GT_SUFFIXERATOR="/opt/home/jong/tools/gt-1.6.2-Linux_x86_64-64bit-complete/bin/./gt"
LTR_FINDER="/opt/home/jong/tools/LTR_FINDER_parallel-1.1/./LTR_FINDER_parallel"
LTR_RETRIEVER="/opt/home/jong/tools/LTR_retriever/LTR_retriever"

# Step 1: Suffixerator
echo "Running gt suffixerator..."
$GT_SUFFIXERATOR suffixerator -db $ASSEMBLY_PATH -indexname $INDEX_NAME -tis -suf -lcp -des -ssp -sds -dna

# Step 2: LTR harvest
echo "Running LTR harvest..."
$GT_SUFFIXERATOR ltrharvest -index $INDEX_NAME -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > ${INDEX_NAME}.harvest.scn

# Step 3: LTR Finder
echo "Running LTR Finder..."
$LTR_FINDER -seq $ASSEMBLY_PATH -threads $THREADS -harvest_out -size 1000000 -time 300

# Step 4: Combine SCN files
echo "Generating combined SCN file..."
cat ${INDEX_NAME}.harvest.scn ${ASSEMBLY_BASENAME}.fasta.finder.combine.scn > ${INDEX_NAME}.rawLTR.scn

# Step 5: LTR retriever
echo "Running LTR Retriever..."
$LTR_RETRIEVER -genome $ASSEMBLY_PATH -inharvest ${INDEX_NAME}.rawLTR.scn -threads $THREADS