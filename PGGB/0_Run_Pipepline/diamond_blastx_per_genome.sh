#!/usr/bin/bash

PROJ_DIR="/opt/home/venice/landrace10_asm5/05_post_pggb/output_all_rice/fasta_annot"
FASTA_DIR="$PROJ_DIR/private_per_genome"
DIAMOND_DIR="$PROJ_DIR/private_per_genome/diamond"
DIAMOND_DB="$PROJ_DIR/diamond/db/uniprot_viridiplantae_canonical.dmnd"
LOG_DIR="$DIAMOND_DIR/logs"

# look for the log directory.
# if there's no log directory, create new log_dir
if [ ! -d "$DIAMOND_DIR" ]; then
	echo "Directory $LOG_DIR does not exist. Creating it now..."
	mkdir -p "$DIAMOND_DIR"
	mkdir -p "$LOG_DIR"
fi

files=($FASTA_DIR/*.fasta)

# Check how many cores we can use
#TOTAL_CORES=$(nproc)
#JOBS=14

# use either of these
#THREADS=$((TOTAL_CORES / JOBS))
THREADS=3

#for fasta in ./private/private.fasta.masked \
#	     ./dispensable/dispensable.fasta.masked \
#	     ./core/core.fasta.masked
for fasta in "${files[@]}"
do
	(
		base=$(basename $fasta .fasta)
		
		echo "Now processing: $fasta...\n"
		#echo "Will write output to $DIAMOND_DIR/${base}_x_uniprot_fmt6.tsv"
		diamond blastx \
		--query "$fasta" \
		--db "$DIAMOND_DB" \
		--evalue 1e-6 \
		--outfmt 6 \
		--max-target-seqs 1 \
		--max-hsps 1 \
		--threads "$THREADS" \
		--out "$DIAMOND_DIR/${base}_x_uniprot_fmt6.tsv"
	) &
done


wait

echo "Diamond Blastx pipeline finished!"
