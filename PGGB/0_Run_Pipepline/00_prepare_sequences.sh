#!/bin/bash

PROJ_DIR="/opt/home/venice/landrace10"  # To do: edit to dynamically change project directory name
ASM_DIR="$PROJ_DIR/assembly_files"	# This directory should contain all the raw assembly files you have (fasta format)
THREADS=8				# Depends on your computer/machine
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")	# 

cd "$PROJ_DIR" || exit 1

# Create all needed directories for running the pipeline
mkdir -p 00_logs 01_prefixed 02_merged 03_mash 04_pggb 05_post_pggb
echo "Created project directories"
echo "----------------------------------"
LOG_DIR="$PROJ_DIR/00_logs"


{
# Copy all the raw data to the working directory
echo "Copying assemblies..."
cp "$ASM_DIR"/* 01_prefixed/
echo "All assembly files now in 01_prefixed/"
echo "----------------------------------"

# Prefix the fasta headers based on PanSpec naming consensus
echo "Prefixing FASTA headers..."
parallel -j $THREADS '
  fastix -p "{/.}#0#" "{}" > "{.}_prefixed.fa"
  ' ::: 01_prefixed/*.fa 01_prefixed/*.fasta
echo "----------------------------------"

# Concatenate all fasta files from all genomes
echo "Concatenating all assemblies..."
cat 01_prefixed/*_prefixed.fa > 02_merged/all_prefixed.fa
echo "Done concatenating."
echo "----------------------------------"

# Sanity check on the sequence headers so far...
echo "Checking a few headers..."
echo "----------------------------------"
grep "^>" 02_merged/all_prefixed.fa | head
echo "Done checking."
echo "----------------------------------"

# Compress the fasta file
echo "Compressing the fasta file..."
bgzip -@ 4 -c 02_merged/all_prefixed.fa > 02_merged/all_prefixed.fa.gz
echo "Fasta file compression done!"
echo "----------------------------------"

# Index the pan-assembly fasta file
echo "Indexing fasta file..."
samtools faidx 02_merged/all_prefixed.fa.gz
echo "Indexing done!"
echo "----------------------------------"

echo "File preparations done!"

} > "$LOG_DIR/pggb_prepseq_${TIMESTAMP}.log" 2>&1 
