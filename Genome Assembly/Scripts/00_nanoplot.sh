# activate preprocessing
# make sure you installed nanoplot to your environment

if [ "$#" -ne 1 ]; then
        echo "Usage: $0 <input_file>"
        exit 1
fi

# extract accession name
input_file=$1
accession_name=$(basename "$input_file" .fastq)

NanoPlot --fastq "$input_file" --outdir "NanoPlot_${accession_name}" --threads 8 --loglength