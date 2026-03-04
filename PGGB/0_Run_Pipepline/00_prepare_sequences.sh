#!/bin/bash

PROJ_DIR="/opt/home/venice/landrace10_asm5" # To do: edit to dynamically change project directory name
ASM_DIR="$PROJ_DIR/assembly_files"
THREADS=8
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

cd "$PROJ_DIR" || exit 1

mkdir -p 00_logs 01_prefixed 02_merged 03_mash 04_pggb 05_post_pggb
echo "Created project directories"
echo "----------------------------------"
LOG_DIR="$PROJ_DIR/00_logs"

{

echo "Copying assemblies..."
cp "$ASM_DIR"/* 01_prefixed/
echo "All assembly files now in 01_prefixed/"
echo "----------------------------------"

echo "Prefixing FASTA headers..."
parallel -j $THREADS '
  fastix -p "{/.}#0#" "{}" > "{.}_prefixed.fa"
  ' ::: 01_prefixed/*.fa 01_prefixed/*.fasta
echo "----------------------------------"

echo "Now preparing one large fasta file."
echo "Concatenating all assemblies..."
cat 01_prefixed/*_prefixed.fa > 02_merged/all_prefixed.fa
echo "Done concatenating."
echo "----------------------------------"

echo "Checking a few headers..."
echo "----------------------------------"
grep "^>" 02_merged/all_prefixed.fa | head
echo "Done checking."
echo "----------------------------------"

echo "Compressing the fasta file..."
bgzip -@ 4 -c 02_merged/all_prefixed.fa > 02_merged/all_prefixed.fa.gz
echo "Fasta file compression done!"
echo "----------------------------------"

echo "Indexing fasta file..."
samtools faidx 02_merged/all_prefixed.fa.gz
echo "Indexing done!"
echo "----------------------------------"

echo "File preparations done!"

} > "$LOG_DIR/pggb_prepseq_${TIMESTAMP}.log" 2>&1 
