# Simulations
In this markdown, I enumerated how to re-create the simulated assemblies that I used to make a comparison between Minigraph, Minigraph-Cactus, and PGGB. I also placed the commands that I used to generate the graphs for these softwares/pipelines. 

Programs:
1. VISOR v1.1.2
2. Bedtools v2.27.1
3. Samtools v1.18
4. Minigraph v0.21 (r606)
5. PGGB v0.7.2
6. vg v1.59.0 "Casatico"
7. ODGI v0.9.0-0-g1895f496
8. Singularity v3.8.6
9. Minigraph-Cactus v2.9.3-gpu

## VISOR
### Steps
1. Fetch BED files from repository and copy the reference genome IRGSP 1.0.
2.Use VISOR to create the simulated assemblies based on the BED files. The code below shows an example of commands to generate the three simulated assemblies including three chromosomes with three variants as dictated in the BED file. 
```
mkdir visor_simulations
# copy bed

VISOR HACk -g /opt/home/jong/oryza/simulation/IRGSP_sample.fasta -b asm1.bed -o asm1.out

VISOR HACk -g /opt/home/jong/oryza/simulation/IRGSP_sample.fasta -b asm2.bed -o asm2.out

VISOR HACk -g /opt/home/jong/oryza/simulation/IRGSP_sample.fasta -b asm3.bed -o asm3.out
```
## Minigraph
### Steps
1. Integrate into Pangenome graph using 16 threads.
```
# integrate
/usr/bin/time -v bash -c '
~/tools/minigraph/./minigraph -cxggs -t16 ../../IRGSP_sample.fasta ../asm1.out/asm1.renamed.fa ../asm2.out/asm2.renamed.fa  ../asm3.out/asm3.renamed.fa  > minigraph_sim_1.gfa
'
```
2. Call the variants.
```
# see bubbles
~/tools/gfatools/./gfatools bubble deletions_1.gfa > deletions_1.tsv

# call bubbles
for asm in ../IRGSP_chr03.fasta deletion_b.fasta deletion_c.fasta deletion_d.fasta deletion_e.fasta; do
basename=$(basename "$asm" .fasta)
  ~/tools/minigraph/./minigraph -cxasm --call deletions_1.gfa $asm > ${basename}.alleles.bed
done
```

## PGGB
### Steps
1. Rename the sample headers to follow panSN-specifications.
```
# rename headers
awk '/^>/{print ">asm1#0#" substr($0, 2); next}1' asm1.fa > asm1.renamed.fa

awk '/^>/{print ">asm2#0#" substr($0, 2); next}1' asm2.fa > asm2.renamed.fa

awk '/^>/{print ">asm3#0#" substr($0, 2); next}1' asm3.fa > asm3.renamed.fa

awk '/^>/{print ">IRGSP#0#" substr($0, 2); next}1' IRGSP_sample.fasta > IRGSP_sample.renamed.fasta
```

2. Merge, zip, and index the file.
```
cat IRGSP_sample.renamed.fasta ../visor_asm_simulations/asm1.out/asm1.renamed.fa ../visor_asm_simulations/asm2.out/asm2.renamed.fa ../visor_asm_simulations/asm3.out/asm3.renamed.fa > asm0123.fa

bgzip -@ 4 asm0123.fa
samtools faidx asm0123.fa.gz
```
3. Calculate divergence and estimate repeats length for -s parameter. 

```
# for -p parameter
mash triangle -s 10000 asm0123.fa.gz > asm0123.10k.mash_triangle.txt
sed 1,1d asm0123.10k.mash_triangle.txt | tr '\t' '\n' | grep chr -v | LC_ALL=C sort -g -k 1nr | uniq | head -n 1

# for -s parameter
awk '{print $1, $5-$4}' repeats01-03.tsv > repeats01-03_diff.tsv
awk 'NR == 1 || $2 > max {max = $2} END {print max}' repeats01-03_diff.tsv
```

4. Run PGGB using 16 threads. 
```
/usr/bin/time -v bash -c ' pggb -i asm0123.fa.gz \ 
    -p 80 -s 25000 \
    -n 4 
    -m 
    -V 'IRGSP' 
    -t 16
    -o asm0123.out
'
```

5. Call the variants using vg deconstruct. Flatten the graph for statistics purposes. 
```
# call the variants
vg deconstruct -P 'IRGSP' asm0123.fa.gz.0f0f8e3.11fba48.ab4226d.smooth.final.gfa --all-snarls --path-traversals --theads 8 --verbose

# flatten the graph
odgi flatten -i asm0123.fa.gz.0f0f8e3.11fba48.ab4226d.smooth.final.gfa -f asm0123.fa.gz.0f0f8e3.11fba48.ab4226d.smooth.final.flattened.fa -t 8
```

## Minigraph-Cactus (MC)
Note that MC requires singularity to run. 
### Steps
1. Create seqfile for MC. 
```
echo "IRGSP /opt/home/jong/oryza/simulation/IRGSP_sample.fasta
asm1 /opt/home/jong/oryza/simulation/visor_asm_simulations/asm1.out/asm1.fa
asm2 /opt/home/jong/oryza/simulation/visor_asm_simulations/asm2.out/asm2.fa
asm3 /opt/home/jong/oryza/simulation/visor_asm_simulations/asm3.out/asm3.fa" > asms.seqfile
```

2. Run MC. Note that MC is a pipeline and inside it, Minigraph uses 4 threads by default. Thus, setting consCores to 4 will set the threads to run Minigraph with 16 threads.
```
mkdir visor_simulation

/usr/bin/time -v bash -c '
singularity exec docker://quay.io/comparative-genomics-toolkit/cactus:v2.9.3-gpu \
cactus-pangenome ./js ./asms.seqfile --outDir ./visor_asms-pg --outName visor_asms --reference IRGSP --filter 2 --haplo --giraffe filter --viz --odgi --chrom-vg clip filter --chrom-og --gbz clip filter --gfa clip full --vcf --vcfReference IRGSP --logFile ./visor_asms.log --workDir ./visor_simulation --consCores 4 --mgMemory 128Gi 
'
```

## Sample Results



| Parameter | Minigraph | Minigraph - Cactus | PGGB |
| --------- | --------- | ---------- | --------- |
| Wall-clock time (h:mm:ss)	| 04:21.7 | 13:27.5	| 18:26.98 |
| System time (seconds) |	16.36 | 163.36 | 1841.66 |
| Memory (Maximum RSS) (Gb) |	6	| 6.416	| 11.223864 |
