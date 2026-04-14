#!/bin/bash

threads=5
engine="ncbi"
repbase_species="rice"
start_time=$(date +%s)

rep_library="/opt/home/venice/landrace10/rerun_output_all_O_mer_p90_noseqpart/fasta_annot/pggb_all_db-families.fa"
dir="/opt/home/venice/landrace10/rerun_output_all_O_mer_p90_noseqpart/fasta_annot/private_fasta"

mkdir -p logs

run_rm () {

	fasta=$1
	base_name=$(basename "$fasta" .fasta)
	output_dir="$dir/$base_name"

	mkdir -p "$output_dir"

	echo "[$(date)] Running RepeatMasker for $fasta"

	RepeatMasker \
	        -e "$engine" \
        	-pa "$threads" \
		-lib "$rep_library" \
		-xsmall \
		-dir "$output_dir" \
		"$dir/$fasta"

	masked_file=$(ls $output_dir/*.masked 2>/dev/null | head -n1)

	if [[ -f "$masked_file" ]]; then
		RepeatMasker \
			-e "$engine" \
			-pa "$threads" \
			-species "$repbase_species" \
			-xsmall \
			-dir "$output_dir" \
			"$masked_file"
	else
		echo "WARNING: masked file not found for $fasta"
	fi

	zcat $output_dir/*.cat.gz | gzip > "$output_dir/${base_name}.full.cat.gz"

	ProcessRepeats \
		-a \
		-species "$repbase_species" \
		"$output_dir/${base_name}.full.cat.gz" \
		> logs/${base_name}.fullmask.log 2>&1
	
	echo "[$(date)] Finished $fasta"
} > logs/run_repeat_masker_parallel_"$start_time".log 

fasta_files=`ls $dir`

# Launch jobs in background
for f in $fasta_files
do
	 run_rm "$dir/$f" & 
done

#run_rm core.fasta &
#run_rm dispensable.fasta &
#run_rm private.fasta &

# Wait until all bg jobs are finished
wait

echo "Pipeline completed"

end_time=$(date +%s)
runtime=$((end_time - start_time))

echo "Total runtime: $runtime seconds"
echo "Total runtime: $((runtime/60)) minutes"
echo "Total runtime: $((runtime/3600)) hours"
