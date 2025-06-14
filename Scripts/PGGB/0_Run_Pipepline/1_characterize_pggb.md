# PGGB Analysis Pipeline
After creating the graph pangenome, the following steps enumerate a workflow for characterizing the pangenome based on sequences encoded on the nodes.

Make sure your current environment as the following programs installed:
1. panacus v.0.2.5
2. odgi v0.9.0-0-g1895f496
3. seqkit v2.6.0
4. RepeatMasker v4.1.7
5. Diamond v2.1.10

## Steps

1. Go to directory of generated pangenomes
```
cd $asm_dir
```
2. Merge the pangenomes per chromosome and characterize by quantifying number of countables and pangenome size using Panacus.
   
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

4. Categorize the countables into the following: Core, dispensable, and private. Categorization is based on the coverage of each asssembly on the countables as determined by Panacus. In this case, the countables are nodes.
```
panacus table --count node --groupby-sample --threads 4 $gfa_asm5 > asm5_graphs.table.tsv

# add sum column
awk 'NR==1 {print $0 "\tSum"; next} {sum=0; for(i=2;i<=NF;i++) sum+=$i; print $0 "\t" sum}' asm5_graphs.table.tsv > asm5_graphs.table.sum.tsv

# quantify node coverage total
seq 0 7 | while read i; do 
awk -v i="$i" '$2 == i' asm5_graphs.table.total.tsv | wc -l; 
done
```
5. Extract node IDs per category. Column number 2 is the total column. A sum of 7 (total number of samples) means presence of all assembles. Therefore, these nodes are considered core. For dispensable, these values become less than 7 or greater than 2. Meanwhile, a value of 1 means private category. Note that these two categories were merged into the "variable" category. The definition of these categories depend on the context of the study.
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

6. Determine quantity of repetitive elements. Edit variables inside `repeat_masker_pggb.sh`
```
bash repeat_masker_pggb.sh
```

After this step, you may choose to visualize the fractions of repeats. See `scripts/0_Visualize_Pangenome_stats_repeats.R `as an example.

7. Blast repeatmasked fasta sequences using
diamond against the UniProt/SwissProt database (viridiplantae only).
```
# create a database
diamond makedb --in ~/db/uniprot_viridiplantae_canonical.fasta -d ~/db/uniprot_viridiplantae_canonical.dmnd

# blast the sequences
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
This step only took less than 5 minutes for a pangenome with 7 samples.

9.  For the annotation of homologous hits, go to https://pantherdb.org for Gene Ontology enrichment. Use PantherDB -> Upload List -> Oryza Gene set -> Enrichment test -> Fisher's exact -> Calculate FDR. Query the blast hits for each category or groups of interest and download the resulting tsv file.

Use a help script called `Scripts/panther_prep.R` to place a sample columnn on the tsv output. For example,
```
Rscript panther_prep.R core.txt cre
```
You should have enrichment results for the following GO:
   - Biological Process
   - Molecular Function
   - Protein Class


Visualize using the following scripts. These scripts are not automated so edit the content as needed. You may also use your own visualization scripts.

| Task | Script |
|------| -------|
| Enrichment visualization | `Scripts/0_Visualize_Categories_enrichment.R` |
| Family and domain analysis | `Scripts/0_GeneFamily_domains_analyses.R`|
|

10. Move to `2_Gain_Loss_analysis.md`
