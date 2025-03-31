#!/bin/bash

# Usage: ./script.sh <genome.fa> <threads> <long_reads.fastq.gz>

# Check if correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <genome.fa> <threads> <long_reads.fastq.gz>"
    exit 1
fi

# Assign input arguments
GENOME=$1
THREADS=$2
LONG_READ=$3

# Activate purge_haplotigs environment (modify as needed)
# source activate purge_haplotigs  # Uncomment if using conda
NAME=$(basename "$GENOME" .fasta)
# Index the assembly
bwa index $GENOME

# Map reads using minimap2
minimap2 -t $THREADS -ax map-pb $GENOME $LONG_READ --secondary=no \
    | samtools sort -m 1G -o "${NAME}.bam" -T tmp.ali

# STEP 1: Generate histogram
purge_haplotigs hist -b "${NAME}.bam" -g $GENOME -t $THREADS