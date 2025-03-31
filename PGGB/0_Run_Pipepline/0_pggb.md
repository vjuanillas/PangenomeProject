# PGGB Pangenome

The following were the steps taken to create the pangenome graph composed of 7 assemblies, including IRGSP 1.0 and O. meridionalis as outgroup.

Programs:
1. PGGB v.0.7.2
2. fastix v.0.1.0
3. mash v.2.3
4. samtools v.1.21
5. bgzip v.1.21

## Steps:

1. Prepare input sequences. Make sure that the name of the fasta files correspond to the sample names.
```
working_dir=pwd

cd $ASM_DIR

sed '/^>/ s/^>\([^ ]*\).*/>out_chr\1/' O_mer.fa > O_mer.renamed.fa

mv cw02.renamed.fasta Or.renamed.fasta
mv nh232.renamed.fasta Osj.renamed.fasta
mv nh236.renamed.fasta Osi.renamed.fasta
mv nh286.renamed.fasta Ob.renamed.fasta
mv nh273.renamed.fasta Og.renamed.fasta

# fastix
for asm in *.renamed.fasta; do
basename=$(basename "$asm" .renamed.fasta);
fastix -p "${basename}#0#" $asm > "${basename}_prefixed.fa"; 
done 

export PATH="/opt/home/jong/.cargo/bin:$PATH"

fastix -p "O_mer#0#" $asm > "O_mer_prefixed.fa" 
done

fastix -p "Osj(IRGSP)#0#" $asm > "IRGSP.prefixed.fa" 
done

# merge all fasta files
cat O_mer_prefixed.fa IRGSP_prefixed.fa Or_prefixed.fa Osj_prefixed.fa Osi_prefixed.fa Ob_prefixed.fa Og_prefixed.fa> all_rice_O_mer.fa

bgzip -@ 4 all_rice_O_mer.fa
samtools faidx all_rice_O_mer.fa.gz


```

2. Calculate mash distance for -p calculation. -p = 100 - max_divergence*100
```
# try two parameters of -s
mash triangle all_rice_O_mer.fa.gz -i > all_rice_O_mer.mash_triangle.default.txt

mash triangle all_rice_O_mer.fa.gz -s 10000 -i > all_rice_O_mer.mash_triangle.10k.txt

sed 1,1d all_rice_O_mer.mash_triangle.default.txt | tr '\t' '\n' | grep chr -v | LC_ALL=C sort -g -k 1nr | uniq | head -n 2

sed 1,1d all_rice_O_mer.mash_triangle.10k.txt | tr '\t' '\n' | grep chr -v | LC_ALL=C sort -g -k 1nr | uniq | head -n 2
```
3. Estimate segment length for pggb based on repeats of assemblies. Repetitive element data from RepeatMasker was used for -s estimation.
```
sed 1,3d asm5.npb.fa.out | awk '{print $6, $7, $7-$6}' | awk '$3 > max {max = $3} END {print max}' | head
```

4. Partition the sequences into communities
```
partition-before-pggb -i all_rice_O_mer.fa.gz \
                      -o output_all_O_mer_p90 \
                      -n 7 \
                      -t 16 \
                      -p 90 \ # replace with 60
                      -s 50k \
                      -V 'O_mer:100000'
```
5. Run PGGB on the 12 communities corresponding to each Oryza chromosomes
```
/usr/bin/time -v bash -c '
seq 0 11 | while read i; do
pggb -i output_all_O_mer_p90/all_rice_O_mer.fa.gz.dac1d73.community.$i.fa \
     -o output_all_O_mer_p90/all_rice_O_mer.fa.gz.dac1d73.community.$i.fa.out \
     -s 50000 -l 250000 -p 90 -c 1 -K 19 -F 0.001 -g 30 \
     -k 23 -f 0 -B 10M \
     -n 7 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ \
     -Y "#" -V O_mer:100000 --threads 16 --poa-threads 16
done
'
```