# PGGB Analysis Pipeline
Programs:
- panacus v.0.2.5
- odgi v0.9.0-0-g1895f496
- seqkit v2.6.0

## Steps

1. Go to directory of generated pangenomes
```
cd $asm_dir
```
2. Merge the pangenomes per chromosome and characterize by quantifying number of countables and pangenome size. 
   
```
ls $(pwd)/all_rice_O_mer.fa.gz.dac1d73.community.{0..11}.fa.out/*smooth.final.og > asm5_og_graphs.txt

# merge
odgi squeeze -f asm5_og_graphs.txt -o asm5_graphs.og -threads=4 

# convert to gfa
for og in *.og; do basename=$(basename "$og" .og); odgi view --idx=${og} --to-gfa > ${basename}.gfa ; done

# get statistics
panacus info asm5_graphs.gfa -S -t 4 > asm5_graphs.info.tsv
```

3. Visualize the pangenome growth, including number of new nodes/bp upon addition of assemblies.
```
# prepare sample names
echo 'Osj(IRGSP) O_mer Or Osj Osi Ob Og' | tr ' ' '\n' > asm5.order.samples.txt

# declare variable
gfa_asm5=~/oryza/asm/pggb_pangenomes/output_all_O_mer_p90/asm5_graphs.gfa

# for bp
RUST_LOG=info panacus ordered-histgrowth -c bp -t4 -l 1,2,7 -S -O asm5.order.samples.txt $gfa_asm5 > asm5.ordered-histgrowth.bp.tsv

# for nodes
RUST_LOG=info panacus ordered-histgrowth -c node -t4 -l 1,2,7 -S -O asm5.order.samples.txt $gfa_asm5 > asm5.ordered-histgrowth.node.tsv

# visualize
panacus-visualize asm5.ordered-histgrowth.bp.tsv > asm5.ordered-histgrowth.bp.pdf

panacus-visualize asm5.ordered-histgrowth.node.tsv > asm5.ordered-histgrowth.node.pdf
```
Sample visualization:
![[Pasted image 20250324131826.png | 400]]

4. Core, dispensable, and private. Calculate coverage of nodes by using Panacus.
```
panacus table --count node --groupby-sample --threads 4 $gfa_asm5 > asm5_graphs.table.tsv

# add sum column
awk 'NR==1 {print $0 "\tSum"; next} {sum=0; for(i=2;i<=NF;i++) sum+=$i; print $0 "\t" sum}' asm5_graphs.table.tsv > asm5_graphs.table.sum.tsv

# quantify node coverage total
seq 0 7 | while read i; do 
awk -v i="$i" '$2 == i' asm5_graphs.table.total.tsv | wc -l; 
done
```

| Number of Assemblies | Number of nodes | Description |
| -------------------- | --------------- | ----------- |
| 1                    | 8609186         | Private     |
| 2                    | 5376727         | Dispensable |
| 3                    | 3466444         | Dispensable |
| 4                    | 3495548         | Dispensable |
| 5                    | 3351433         | Dispensable |
| 6                    | 6210558         | Dispensable |
| 7                    | 3850484         | Core        |
| TOTAL                | 34 360 380      |             |
| DISP TOTAL           | 21 900 710      |             |
5. Extract node IDs per category
Tool: seqkit grep
```
awk '$2 == 7' asm5_graphs.table.total.tsv | cut -f1 > asm5_graphs.core.txt

awk '$2 == 1' asm5_graphs.table.total.tsv | cut -f1 > asm5_graphs.private.txt

awk '$2 > 1 && $2 < 7' asm5_graphs.table.total.tsv | cut -f1 > asm5_graphs.dispensable.txt

# declare variable
PAN_FASTA=/opt/home/jong/oryza/asm/pggb_pangenomes/output_all_O_mer_p90/asm5_graphs.fa
DIR_ANNOT=/opt/home/jong/oryza/asm/pggb_pangenomes/output_all_O_mer_p90/fasta_annot

for text in asm5_graphs*.txt; do
basename=$(basename "$text" .txt)
seqkit grep -f ${text} "${PAN_FASTA}" > "${DIR_ANNOT}/${basename}.fasta"
done

# do the same for private 
# adjust to node columns

for col in {2..8}; do
    outfile=$(awk -v c="$col" 'NR==1 {print $c}' ../asm5_graphs.table.sum.tsv)  # Get the header name
    awk -v c="$col" '$9 == 1 && $c == 1 {print $1}' ../asm5_graphs.table.sum.tsv > "asm5_graphs.${outfile}.tsv"
done
```

6. Determine quantity of repetitive elements. Edit variables inside repeat_masker_pggb.sh
```
bash repeat_masker_pggb.sh
```

Visualize

7. Blast repeatmasked fasta sequences
diamond version 2.1.10
```
diamond makedb --in ~/db/uniprot_viridiplantae_canonical.fasta -d ~/db/uniprot_viridiplantae_canonical.dmnd

for fasta in asm5_graphs.core.fasta.masked asm5_graphs.dispensable.fasta.masked asm5_graphs.private.fasta.masked; do
    basename=$(basename "$fasta" .fasta.masked)
    diamond blastx \
        --query "$fasta" \
        --db ~/db/uniprot_viridiplantae_canonical.dmnd \
        --evalue 1e-6 \
        --outfmt 6 \
        --max-target-seqs 1 \
        --max-hsps 1 \
        --threads 16 \
        --out "${basename}_x_uniprot_fmt6.tsv"
done


```

8. For each hits in dispensable and private, map node IDs and extract blast hits
```
# take the column with query node number from blast tsv

# DECLARE
DIR=/opt/home/jong/oryza/asm/pggb_pangenomes/output_all_O_mer_p90/fasta_annot/

for tsv in ${DIR}/asm5_graphs.private/*.tsv; do
    basename=$(basename "$tsv" .tsv)
    awk '{split($2, a, "|"); print $1, a[3]}' "$tsv" > "${basename}_clean.tsv"
done

# given node numbers, take the matching node on Column 1 


GENE_LIST=/opt/home/jong/oryza/asm/pggb_pangenomes/output_all_O_mer_p90/fasta_annot/asm5_graphs.private/asm5_graphs.private_x_uniprot_fmt6_diamond_clean.tsv

PRIV_GRPS=/opt/home/jong/oryza/asm/pggb_pangenomes/output_all_O_mer_p90/panacus/private_groups

OUT=private_hits
mkdir private_hits

for tsv in $PRIV_GRPS/*.tsv; do
basename=$(basename "$tsv" .tsv)
awk 'NR==FNR {nodes[$1]; next} $1 in nodes' $tsv $GENE_LIST > ${basename}.matched.tsv
awk '{print $2}' ${basename}.matched.tsv > ${OUT}/${basename}.hits.tsv
done
```

9.  Use PantherDB -> Upload List -> Oryza Gene set -> Enrichment test -> Fisher's exact -> Calculate FDR. Manually place a new column called "Sample" using any TSV viewer (e.g. Excel)
   - Biological Process
   - Molecular Function
   - Protein Class
Visualize using the script
`~/Desktop/LAB_Files/pggb/editing_scripts/0_Visualize_Categories_enrichment.R`