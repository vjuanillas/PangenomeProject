# Genome assembly of Taiwanese accessions
This is a summary of how to assemble ONT simplex reads. 

## Tools:
1. Dorado v.0.8.0
2. Samtools v.1.20
3. Minimap2 v.2.28-r1209
4. Nanoplot v.1.42.0
5. Chopper v.0.8.0
6. seqkit v.2.6.0
7. Flye v.2.9.5-b1801
8. tgsgapcloser v.1.2.1
9. Medaka v.1.11.3
10. BUSCO v.5.7.1
11. To compute for LAI, LTR_retriever (v2.9.8) package is needed.


## Steps:

1. Perform basecalling on preferred models. Make sure that the models and version of Dorado are consistent all throughout with your work. You may opt to add --emit-fasta if BAM files are not necessary. 
```
/opt/home/jong/tools/dorado-0.8.0-linux-x64/bin/dorado basecaller sup --min-qscore 10 /opt/data/jong/minION/AS132/pod5_2/pod5_0423 > AS132_2_simplex.bam

samtools fastq AS132_2_simplex.bam > AS132_2_simplex.fastq

# there were two runs for this accession. Hence, it was combined simply using the cat command as follows
cat AS132_simplex_1.fastq AS132_2_simplex.fastq > AS132_simplex_12.fastq
```

Optional: You can calculate the coverage by aligning the reads to IRGSP 1.0 reference genome and make sequences visualization using nanoplot. All of these processes can be automated using two scripts deposited in `/Scripts/`
```

# nanoplot combined
bash ../../scripts/00_nanoplot.sh /opt/data/jong/minION/AS132/00_Raw_Reads/AS132_simplex_12.fastq

# coverage
bash ../../scripts/00_coverage.sh /opt/data/jong/minION/AS132/00_Raw_Reads/AS132_simplex_12.fastq
```

2. The raw sequences can be preprocessed using chopper, which combines both filtering shorter reads and low quality sequences. 
```
# preprocess
bash ../../scripts/00_preprocessing.sh /opt/data/jong/minION/AS132/00_Raw_Reads/AS132_simplex_12.fastq
```

3. The sequences can be corrected using the HERRO algorithm. This can be done through Dorado. Note that HERRO will filter reads with less 10kb in length and less than 10 q score. Sometimes, there are reads that are just 1 bp in length. These can be filtered out using seqkit.
```
/opt/home/jong/tools/dorado-0.8.0-linux-x64/bin/dorado correct -t 8 /opt/data/jong/minION/AS132/01_Preprocessing/AS132_simplex_12.preprocessed.fastq > AS132_simplex.corrected.fasta

# filter out short reads
cat AS132_simplex.corrected.fasta | seqkit seq -m 1000 > AS132_simplex.corrected.filtered.fasta
```

4. You can then assemble the corrected sequences using Flye. Note that scaffolding is enabled here. You may opt to disable this if not necessary.
```
flye --nano-corr /opt/data/jong/minION/AS132/01_Preprocessing/AS132_simplex.corrected.filtered.fasta --genome-size 400m -o 1_AS132_flye --threads 8 --scaffold
```

5. If scaffolding was done, gaps will be represented as strings of "N" or "n". Here, you can close the gaps with the raw sequencing reads using tgsgapcloser.
```
tgsgapcloser --scaff /opt/data/jong/minION/AS132/02_Assembly/1_AS132_flye/assembly.fasta --reads /opt/data/jong/minION/AS132/01_Preprocessing/AS132_simplex.corrected.filtered.fasta --output 02_AS132_gclosed --ne > 2_AS132_gclosed.log 2>2_AS132_gclosed.err
```

6. You can polish the assembly using medaka. Note that medaka requires the correct model used to basecall the reads from the pod5.
```
medaka_consensus -i /opt/data/jong/minION/AS132/00_Raw_Reads/AS132_simplex_12.fastq -d /opt/data/jong/minION/AS132/02_Assembly/2_gclosed/02_AS132_gclosed.scaff_seqs -o 3_polishing -t 8 -m 'r1041_e82_400bps_sup_v5.0.0'
```

7. You can calculate BUSCO and lai for quality assessment. 

```
# BUSCO
busco -c 4 -m genome -i /opt/data/jong/minION/AS132/02_Assembly/3_polishing/consensus.fasta -o busco_AS132_polished --lineage_dataset embryophyta_odb10

# lai
mkdir lai
bash /opt/data/jong/minION/scripts/01_lai.sh /opt/data/jong/minION/AS132/02_Assembly/3_polishing/consensus.fasta AS132 8
```



