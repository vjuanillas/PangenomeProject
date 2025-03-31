#!/bin/bash/

# activate conda pangenome_env
# make sure to have minimap2 and samtools installed
# check if correct number of arguments
if [ "$#" -ne 1 ]; then
        echo "Usage: $0 <input_file>"
        exit 1
fi

# extract accession name
input_file=$1
accession_name=$(basename "$input_file" .fastq)

# align to reference genome
minimap2 -x map-ont -a ~/oryza/reference_genome/IRGSP-1.0_genome.fasta "$input_file" | samtools sort -o "${accession_name}_ref_mapont.sorted.bam" -@ write-index

samtools depth -aa "${accession_name}_ref_mapont.sorted.bam"  > "${accession_name}_depth.txt"
awk '{sum+=$3} END {print "Mean depth = ", sum/NR}' "${accession_name}_depth.txt"