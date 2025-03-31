#!/bin/bash
GENOMEFA="1_pilon_assembly.fa"
OUTDIR=pilon_out_2
PILONJAR=/opt/home/jong/miniconda3/envs/pangenome_env/share/java/pilon-1.18-0/pilon.jar
JAVA_TOOL_OPTIONS="-Xmx64G -Xss2560k" # set maximum heap size to something reasonable
####### Split the genome
# split the genome into single chromsomes 
bioawk -c fastx '{print $name}' $GENOMEFA > nms
mkdir -p split
samtools faidx $GENOMEFA
for CHR in $(cat nms);do 
        samtools faidx $GENOMEFA $CHR > split/$CHR\.fa; 
done
# this produced split/ directory with individual fastas for each chromosome
####### Run pilon
mkdir -p $OUTDIR
# increase maximum java heap size, as the application crashes otherwise.
# Other option is to split the genome into chunks corresponding to single chromsomes
ls split/*fa > toprocess

for CHRFILE in $(cat toprocess); do
        echo $CHRFILE
        CHRNAME=$(basename $CHRFILE | cut -f 1 -d '.')
        CHROUTDIR=$OUTDIR/$CHRNAME

java $JAVA_TOOL_OPTIONS -jar "$PILONJAR" --nostrays --changes --vcf --verbose --threads 8 --genome $CHRFILE --output $CHRNAME --outdir $CHROUTDIR  --frags bwa_mapping_2_sorted.bam >> pilon_2.logs 2>&1
done