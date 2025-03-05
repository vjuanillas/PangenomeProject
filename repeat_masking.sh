#!/bin/bash

# Define the repeat library
rep_library="/opt/home/jong/oryza/asm/minigraph_pangenomes/00_mash_distance/repeatmodeler/library/asm5-families.fa"
# indels
# Define the directory containing the FASTA files
dir="/opt/home/jong/oryza/asm/pggb_pangenomes/output_all_O_mer_p90/fasta_annot/indel"
# List of FASTA files
fasta_files=(
    "gain_indel_africa.fasta" "gain_indel_japonicWild.fasta" "gain_indel_Or.fasta"
    "gain_indel_Osj_ref.fasta" "loss_indel_japonica.fasta" "loss_indel_Og.fasta"
    "loss_indel_Osj.fasta" "gain_indel_asia.fasta" "gain_indel_Ob.fasta"
    "gain_indel_Osi.fasta" "loss_indel_africa.fasta" "loss_indel_japonicWild.fasta"
    "loss_indel_Or.fasta" "loss_indel_Osj_ref.fasta" "gain_indel_japonica.fasta"
    "gain_indel_Og.fasta" "gain_indel_Osj.fasta" "loss_indel_asia.fasta"
    "loss_indel_Ob.fasta" "loss_indel_Osi.fasta"
)

# Loop through each FASTA file
for fasta in "${fasta_files[@]}"; do
    # Get the basename without the .fasta extension
    base_name=$(basename "$fasta" .fasta)

    # Define the output directory for this file
    output_dir="$dir/$base_name"

    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Run RepeatMasker
    RepeatMasker -e ncbi -pa 16 -lib "$rep_library" -xsmall -dir "$output_dir" "$dir/$fasta"
    echo "First RepeatMasker run completed for $fasta using $rep_library"

    # Find the masked file generated
    masked_file="$output_dir/$(basename "$fasta").masked"

    # Check if the masked file exists before proceeding
    if [[ -f "$masked_file" ]]; then
        echo "Re-masking $masked_file with Repbase library (-species rice)..."

        # Run RepeatMasker again using Repbase (-species rice)
        RepeatMasker -e ncbi -pa 16 -species rice -xsmall -dir "$output_dir" "$masked_file"

        echo "Second RepeatMasker run completed for $masked_file using Repbase (-species rice)"
    else
        echo "Warning: Masked file $masked_file not found, skipping second RepeatMasker run."
    fi
done



# sv
# Define the repeat library
rep_library="/opt/home/jong/oryza/asm/minigraph_pangenomes/00_mash_distance/repeatmodeler/library/asm5-families.fa"
# indels
# Define the directory containing the FASTA files
dir="/opt/home/jong/oryza/asm/pggb_pangenomes/output_all_O_mer_p90/fasta_annot/sv"
# List of FASTA files
fasta_files=(
    "gain_sv_africa.fasta" "gain_sv_japonicWild.fasta" "gain_sv_Or.fasta"
    "gain_sv_Osj_ref.fasta" "loss_sv_japonica.fasta" "loss_sv_Og.fasta"
    "loss_sv_Osj.fasta" "gain_sv_asia.fasta" "gain_sv_Ob.fasta"
    "gain_sv_Osi.fasta" "loss_sv_africa.fasta" "loss_sv_japonicWild.fasta"
    "loss_sv_Or.fasta" "loss_sv_Osj_ref.fasta" "gain_sv_japonica.fasta"
    "gain_sv_Og.fasta" "gain_sv_Osj.fasta" "loss_sv_asia.fasta"
    "loss_sv_Ob.fasta" "loss_sv_Osi.fasta"
)

# Loop through each FASTA file
for fasta in "${fasta_files[@]}"; do
    # Get the basename without the .fasta extension
    base_name=$(basename "$fasta" .fasta)

    # Define the output directory for this file
    output_dir="$dir/$base_name"

    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Run RepeatMasker
    RepeatMasker -e ncbi -pa 16 -lib "$rep_library" -xsmall -dir "$output_dir" "$dir/$fasta"
    echo "First RepeatMasker run completed for $fasta using $rep_library"

    # Find the masked file generated
    masked_file="$output_dir/$(basename "$fasta").masked"

    # Check if the masked file exists before proceeding
    if [[ -f "$masked_file" ]]; then
        echo "Re-masking $masked_file with Repbase library (-species rice)..."

        # Run RepeatMasker again using Repbase (-species rice)
        RepeatMasker -e ncbi -pa 16 -species rice -xsmall -dir "$output_dir" "$masked_file"

        echo "Second RepeatMasker run completed for $masked_file using Repbase (-species rice)"
    else
        echo "Warning: Masked file $masked_file not found, skipping second RepeatMasker run."
    fi
done






