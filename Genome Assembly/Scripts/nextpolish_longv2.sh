#!/bin/bash
logfile="nextpolish_log_$(date +%Y-%m-%d_%H-%M-%S).txt"

echo "Start: $(date "+%Y-%m-%d %H:%M:%S")" | tee -a $logfile

round=2
threads=8
read=../SRR12451701.fastq
read_type=ont
mapping_option=map-ont
input=../nextdenovo_rice/03.ctg_graph/nd.asm.fasta

for i in $(seq 1 $round); do
    echo "Iteration $i Start: $(date "+%Y-%m-%d %H:%M:%S")" | tee -a $logfile

    /usr/bin/time -o polish_time_output_$i.txt -a -f "\nIteration $i: Elapsed Time (hh:mm:ss): %E\n" \
    minimap2 -ax ${mapping_option} -t $threads $input $read | \
    samtools sort - -m 2g --threads $threads -o lgs.sort.bam && \
    samtools index lgs.sort.bam && \
    ls "$(pwd)"/lgs.sort.bam > lgs.sort.bam.fofn && \
    python /opt/home/jong/tools/NextPolish/lib/nextpolish2.py -g $input -l lgs.sort.bam.fofn -r $read_type \
    -p $threads -sp -o genome.nextpolish.fa > polish_$i_log.txt 2>&1 && \

    if [ $i -ne $round ]; then
        mv genome.nextpolish.fa genome.nextpolishtmp.fa
        input=genome.nextpolishtmp.fa
    fi
    echo "Iteration $i End: $(date "+%Y-%m-%d %H:%M:%S")" | tee -a $logfile
done

echo "End: $(date "+%Y-%m-%d %H:%M:%S")" | tee -a $logfile

