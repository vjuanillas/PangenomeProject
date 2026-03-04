 $(pwd)/all_rice_O_mer.fa.gz.dac1d73.community.{0..11}.fa.out/*smooth.final.og > asm5_og_graphs.txt

 # merge
 odgi squeeze -f asm5_og_graphs.txt -o asm5_graphs.og --threads=4 

 # convert to gfa
 for og in *.og; do basename=$(basename "$og" .og); odgi view --idx=${og} --to-gfa > ${basename}.gfa ; done

 # get graph statistics
 panacus info asm5_graphs.gfa -S -t 4 > asm5_graphs.info.tsv



 # prepare sample names
 echo 'Osj(IRGSP) O_mer Or Osj Osi Ob Og' | tr ' ' '\n' > asm5.order.samples.txt

 # declare variable
 gfa_asm5=~/output_all_O_mer_p90/asm5_graphs.gfa

 # for bp
 RUST_LOG=info panacus ordered-histgrowth -c bp -t4 -l 1,2,7 -S -O asm5.order.samples.txt $gfa_asm5 > asm5.ordered-histgrowth.bp.tsv

 # for nodes
 RUST_LOG=info panacus ordered-histgrowth -c node -t4 -l 1,2,7 -S -O asm5.order.samples.txt $gfa_asm5 > asm5.ordered-histgrowth.node.tsv

 # visualize
 panacus-visualize asm5.ordered-histgrowth.bp.tsv > asm5.ordered-histgrowth.bp.pdf

 panacus-visualize asm5.ordered-histgrowth.node.tsv > asm5.ordered-histgrowth.node.pdf

 panacus table --count node --groupby-sample --threads 4 $gfa_asm5 > asm5_graphs.table.tsv

 # add sum column
# awk 'NR==1 {print $0 "\tSum"; next} {sum=0; for(i=2;i<=NF;i++) sum+=$i; print $0 "\t" sum}' asm5_graphs.table.tsv > asm5_graphs.table.sum.tsv

# get the last column number (Sum)
SUM_COL=$(awk 'NR==1{print NF}' asm5_graphs.table.sum.tsv)

# counts for coverage 0..7
for i in $(seq 0 7); do
	awk -v i="$i" -v sum="$SUM_COL" 'NR>1 && $sum==i' asm5_graphs.table.sum.tsv | wc -l
done


 # quantify node coverage total
# seq 0 7 | while read i; do 
# awk -v i="$i" '$2 == i' asm5_graphs.table.total.tsv | wc -l; 
#done

awk '$2 == 7' asm5_graphs.table.total.tsv | cut -f1 > asm5_graphs.core.txt

awk '$2 == 1' asm5_graphs.table.total.tsv | cut -f1 > asm5_graphs.private.txt

awk '$2 > 1 && $2 < 7' asm5_graphs.table.total.tsv | cut -f1 > asm5_graphs.dispensable.txt

# use GFA tools to obtain linearized fasta from graph
gfatools linearize asm5_graphs.gfa > asm5_graphs.fa

# declare variable
PAN_FASTA=~/output_all_O_mer_p90/asm5_graphs.fa
DIR_ANNOT=~/output_all_O_mer_p60/fasta_annot


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


