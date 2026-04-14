#!/bin/bash

PROJ_DIR="/opt/home/venice/landrace10_asm5/05_post_pggb/output_all_rice/fasta_annot/diamond"
LOG_DIR="$PROJ_DIR/logs"
mkdir -p "$LOG_DIR"


echo "Extracting UniProt IDs..."

for group in core dispensable private
do
	infile="${group}_x_uniprot_fmt6.tsv"
	outfile="${group}_uniprot_ids.txt"

	# get the second column, get the middle token which corresponds to the UniProtKB accession id, sort and removethe duplicates 
	cut -f2 "$infile" | \
	cut -d'|' -f2 | \
	sort -u > "$outfile"

	echo "$group: $(wc -l < "$outfile") unique UniProt IDs"
done

# Split the IDs into chunks of 1000 to avoid overwhelming PantherDB 
chunk_file() {
	split -l 1000 -d --additional-suffix=.txt "$1" "$2"
}

#DATABASES=("GO:0008150" "GO:0003674" "PANTHER_PC")

echo "Running PANTHER overrepresentation analysis..."

for group in core dispensable private
do
	LOG_FILE="${LOG_DIR}/${group}_panther.log"
	OUT_DIR="$PROJ_DIR/${group}_chunks"
	mkdir -p $OUT_DIR
	
	IDS=$(paste -sd, ${group}_uniprot_ids.txt)

	echo "[$(date)] Processing $group..." | tee -a "$LOG_FILE"

	chunk_file "${group}_uniprot_ids.txt" "${group}_chunk_"

	for chunk in ${group}_chunk_*.txt
	do
		IDS=$(tr -d '\r' < "$chunk" | paste -sd,)

		base="${chunk%.txt}"
		for dataset in "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP" "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF" "ANNOT_TYPE_ID_PANTHER_PC"
		do
			case $dataset in
				*BP*) suffix="bp" ;;
				*MF*) suffix="mf" ;;
				*PC*) suffix="pc" ;;
			esac

			OUT_FILE="$OUT_DIR/${base}_${suffix}.json"

			echo "Running $dataset on $chunk" | tee -a "$LOG_FILE"

	#curl -X POST "https://pantherdb.org/services/oai/pantherdb/enrich/overrep" \
	#	-H "Content-Type: application/x-www-form-urlencoded" \
	#	-d "geneInputList=${IDS}" \
	#	-d "organism=39947" \
	#	-d "annotDataSet=ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP" \                #Biological Process
	#	-d "enrichmentTestType=FISHER" \
	#	-d "correction=FDR" \
	#	-w "\nHTTP_STATUS:%{http_code}\nTIME_TOTAL:%{time_total}\n" \
	#	-o ${group}_panther_enrichment_bp.json \
	#	2>> "$LOG_FILE" | tee -a "$LOG_FILE"
			
			curl -sSL -X POST "https://pantherdb.org/services/oai/pantherdb/enrich/overrep" \
				-H "Content-Type: application/x-www-form-urlencoded" \
				--data-urlencode "geneInputList=${IDS}" \
				--data-urlencode "organism=39947" \
				--data-urlencode "annotDataSet=${dataset}" \
				--data-urlencode "enrichmentTestType=FISHER" \
				--data-urlencode "correction=FDR" \
				-w "\nHTTP_STATUS:%{http_code}\nTIME_TOTAL:%{time_total}\n" \
				-o "$OUT_FILE" \
				2>> "$LOG_FILE" | tee -a "$LOG_FILE"

			echo "Saved output to $OUT_FILE" | tee -a "$LOG_FILE"
			echo "----------------------------------------" | tee -a "$LOG_FILE"

			sleep 1
		done
	done
done

echo "Finished getting overrep results from PantherDB"
