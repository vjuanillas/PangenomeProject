#!/bin/bash

# ============================= #
# RepeatMasker Pipeline Script  #
# ============================= #

# ------------ Configurable Variables ------------

# Program paths (customize if needed) 
# hashed if with conda; unhash if no conda
# repeatmasker="RepeatMasker"
# processrepeats="ProcessRepeats"

# Number of threads
threads=16

# RepeatMasker options
engine="ncbi"
repbase_species="rice"

# Custom Repeat Library
rep_library="/opt/home/jong/oryza/asm/pggb_pangenomes/output_all_O_mer_p90/RM_lib/pggb_asm5-families.fa"

# Input & Output Directory
dir="/opt/home/jong/oryza/asm/pggb_pangenomes/output_all_O_mer_p90/fasta_annot"

# FASTA files to process
fasta_files=("asm5_graphs.core.fasta" "asm5_graphs.dispensable.fasta" "asm5_graphs.private.fasta")

mkdir -p logs
# ------------ Pipeline ------------

for fasta in "${fasta_files[@]}"; do

    base_name=$(basename "$fasta" .fasta)
    output_dir="$dir/$base_name"

    mkdir -p "$output_dir"

    echo "[$(date)] Running first RepeatMasker for $fasta..."
    RepeatMasker -e "$engine" -pa "$threads" -lib "$rep_library" -xsmall -dir "$output_dir" "$dir/$fasta"
    echo "[$(date)] First RepeatMasker completed for $fasta"

    masked_file="$output_dir/${fasta}.masked"

    if [[ -f "$masked_file" ]]; then
        echo "[$(date)] Re-masking $masked_file using Repbase (-species $repbase_species)..."
        RepeatMasker -e "$engine" -pa "$threads" -species "$repbase_species" -xsmall -dir "$output_dir" "$masked_file"
        echo "[$(date)] Second RepeatMasker completed for $masked_file"
    else
        echo "[$(date)] WARNING: Masked file $masked_file not found, skipping second RepeatMasker run."
    fi
    # ------------ Combine & Post-process ------------

    echo "[$(date)] Concatenating RepeatMasker outputs..."

    cat $output_dir/*cat.* > "$output_dir/${base_name}.full.cat.gz"

    echo "[$(date)] Running ProcessRepeats..."

    ProcessRepeats -a -species "$repbase_species" "$output_dir/${base_name}.full.cat.gz" 2>&1 | tee logs/${base_name}.fullmask.log
done


echo "[$(date)] Pipeline completed!"

# =====================