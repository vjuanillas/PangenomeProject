#!/bin/bash
logfile="sgs_nextpolish_log_$(date +%Y-%m-%d_%H-%M-%S).txt"
echo "Start: $(date "+%Y-%m-%d %H:%M:%S")" | tee -a "$logfile"

#Set input and parameters
round=2
threads=16
read1=/opt/home/jong/oryza/SRR13453606_1.fastq.gz
read2=/opt/home/jong/oryza/SRR13453606_2.fastq.gz
input=lgs_genome.nextpolish.fa

for i in $(seq 1 $round); do
   echo "Iteration $i Start: $(date "+%Y-%m-%d %H:%M:%S")" | tee -a "$logfile"
   /usr/bin/time -o polish_time_output_$i.txt -a -f "\nIteration $i: Elapsed Time (hh:mm:ss): %E\n" \
#step 1:
   #index the genome file and do alignment
   bwa index $input;
   bwa mem -t $threads $input $read1 $read2 |samtools view --threads 3 -F 0x4 -b -|samtools fixmate -m --threads 3  - -|samtools sort -m 2g --threads 5 -|samtools markdup --threads 5 -r - sgs.sort.bam
   #index bam and genome files
   samtools index -@ $threads sgs.sort.bam;
   samtools faidx $input;
   #polish genome file
   python /opt/home/jong/tools/NextPolish/lib/nextpolish1.py -g $input -t 1 -p $threads -s sgs.sort.bam > genome.polishtemp.fa;
   input=genome.polishtemp.fa;
#step2:
   #index genome file and do alignment
   bwa index $input;
   bwa mem -t $threads $input $read1 $read2|samtools view --threads 3 -F 0x4 -b -|samtools fixmate -m --threads 3  - -|samtools sort -m 2g --threads 5 -|samtools markdup --threads 5 -r - sgs.sort.bam
   #index bam and genome files
   samtools index -@ $threads sgs.sort.bam;
   samtools faidx $input;
   #polish genome file
   python /opt/home/jong/tools/NextPolish/lib/nextpolish1.py -g $input -t 2 -p $threads -s sgs.sort.bam > genome.nextpolish.fa;
   input=genome.nextpolish.fa;
   
   echo "Iteration $i End: $(date "+%Y-%m-%d %H:%M:%S")" | tee -a "$logfile"
done;
#Finally polished genome file: genome.nextpolish.fa
echo "End: $(date "+%Y-%m-%d %H:%M:%S")" | tee -a "$logfile"