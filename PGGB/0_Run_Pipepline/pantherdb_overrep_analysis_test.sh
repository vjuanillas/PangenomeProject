#!/bin/bash

# This script uses teh PantherDB API service for overrepresentation/enrichment analysis.
# This reads directly the results from Diamond BLASTx 

PROJ_DIR="/opt/home/venice/landrace10_asm5/05_post_pggb/output_all_rice/fasta_annot/diamond"
LOG_DIR="$PROJ_DIR/logs"
mkdir -p "$LOG_DIR"


# Extract the UniProt IDs
echo "Extracting UniProt IDs..."

# Loop into different categories
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

# Run Overrep Analysis using Fisher Test
echo "Running PANTHER overrepresentation analysis..."


# Get the results per category
for group in core dispensable private
do
	LOG_FILE="${LOG_DIR}/${group}_panther.log"
	OUT_FILE="${group}_panther_enrichment_bp.json"
	
	IDS=$(paste -sd, ${group}_uniprot_ids.txt)

	echo "[$(date)] Processing $group..." | tee -a "$LOG_FILE"

	# For more information regarding the API:
	# https://pantherdb.org/services/tryItOut.jsp?url=%2Fservices%2Fapi%2Fpanther

	curl -X POST "https://pantherdb.org/services/oai/pantherdb/enrich/overrep" \
		-H "Content-Type: application/x-www-form-urlencoded" \
		-d "geneInputList=${IDS}" \
		-d "organism=39947" \			# Oryza sativa
		-d "annotDataSet=GO:0008150" \          # Biological Process (description": "Gene Ontology Biological Process annotations including both manually curated and electronic annotations.)
		#-d "annotDataSet=ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP" \                #Biological Process (GO Slim)
		-d "enrichmentTestType=FISHER" \
		-d "correction=FDR" \
		-w "\nHTTP_STATUS:%{http_code}\nTIME_TOTAL:%{time_total}\n" \
		-o ${group}_panther_enrichment_bp.json \
		2>> "$LOG_FILE" | tee -a "$LOG_FILE"
			

	echo "Saved output to $OUT_FILE" | tee -a "$LOG_FILE"
	echo "----------------------------------------" | tee -a "$LOG_FILE"

	sleep 10
done

echo "Finished getting overrep results from PantherDB"
