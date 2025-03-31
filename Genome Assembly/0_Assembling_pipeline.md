# =============== #
# Genome Assembly #
# =============== #

The following were the steps taken to create the genome assemblies used in this study. 

Programs:
1. NextDenovo v2.5.2
2. NextPolish v1.4.1
3. Purge_haplotigs v1.1.2
4. FCS v0.5.0
5. BUSCO v5.7.1
6. Seqkit v.
7. Mummer NUCmer v3.1
8. gt (LAI) v1.6.2
9. LTR_FINDER v1.0.7
10. LTR_retriever v2.9.8
11. InterProScan version 5.55-88.0

Steps:
1. Preprocess
2. Place in input.fofn the pre-processed sequences to be assembled.
```
ls <FASTA> > input.fofn
```

3. Create the config file for NextDenovo program. Refer to nextdenovo.cfg for the specifications of parameters. Run the program. 
```
logfile="nextden_cor_log_$(date +%Y-%m-%d_%H-%M-%S).txt"

echo "Start: $(date "+%Y-%m-%d %H:%M:%S")" | tee -a "$logfile"

/usr/bin/time -o time_output.txt nextDenovo run.cfg 2>&1 | tee -a "$logfile" 

echo "End: $(date "+%Y-%m-%d %H:%M:%S")" | tee -a "$logfile"
```
The program will output the raw assembly file to be polished using both long and short reads (if available). 

4. Prepare the input file for NextPolish, create the config file, and run the program. 
```
# long reads
ls /opt/home/jong/oryza/ir64_trimmed.fastq > lgs.fofn

# short reads
ls ~/oryza/SRR13453606_1.fastq.gz ~/oryza/SRR13453606_2.fastq.gz > sgs.fofn

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
```
5. 
