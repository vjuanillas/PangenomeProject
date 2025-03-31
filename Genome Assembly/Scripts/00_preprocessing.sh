#!/bin/bash/

# check if correct number of arguments were provided
# make sure you installed porechop_abi and chopped
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# extract accession name from input file
input_file=$1
accession_name=$(basename "$input_file" .fastq)

#porechop
porechop_abi --ab_initio -i "$input_file" -o "${accession_name}.chopped.fastq" -t 8 > porechop.txt 2>&1

# filter short and low quality reads
chopper -q 10 -l 1000 --threads 8 -i "${accession_name}.chopped.fastq" > "${accession_name}.preprocessed.fastq"