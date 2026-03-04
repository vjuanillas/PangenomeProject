#!/bin/bash

PROJ_DIR="$HOME/landrace10"
PGGB_DIR="$PROJ_DIR/04_pggb"
POST_PGGB_DIR="$PROJ_DIR/05_post_pggb"
LOG_DIR="$PROJ_DIR/00_logs"

PART_DIR="$PGGB_DIR/rerun_output_all_rice"
POSTPGGB_OUT_DIR="$POST_PGGB_DIR/rerun_output_all_rice"
mkdir $POSTPGGB_OUT_DIR
THREADS=10
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

{
echo "Merging all PGGB generated OG files..."
find $PART_DIR/ -type f -name "*smooth.final.og" > "$POST_PGGB_DIR/asm_og_graphs.txt"

OG_COUNT=$(wc -l < "$POST_PGGB_DIR/asm_og_graphs.txt")
echo "Found $OG_COUNT OG graphs."

# ------------- Merge all OG files ------------- #
echo "Merging all OG file from each communities..."
odgi squeeze -f "$POST_PGGB_DIR/asm_og_graphs.txt" -o "$POST_PGGB_DIR/merged_graph.og" -t "$THREADS"
echo "Done merging all OG files!"
echo "------------------------------------------------"

# ---------------------------------------------- #
echo "Converting OG file to GFA..."
odgi view -i "$POST_PGGB_DIR/merged_graph.og" --to-gfa > "$POST_PGGB_DIR/merged_graph.gfa"
echo "Done converting OF file to GFA!"
echo "------------------------------------------------"

# ---------------------------------------------- #
echo "Generate linear fasta file from graph..."
gfatools stat "$POST_PGGB_DIR/merged_graph.gfa" > "$POST_PGGB_DIR/merged_graph_gfatool_stat.tsv"
gfatools gfa2fa -l 60 "$POST_PGGB_DIR/merged_graph.gfa" > "$POST_PGGB_DIR/merged_graph_linear_seq.fasta"
echo "Done! Linear fasta sequence written in: $POST_PGGB_DIR/merged_graph_linear_seq.fasta"
echo "------------------------------------------------"

# ---------------------------------------------- #
echo "Generating Graph statistics from merged GFA..."
panacus info "$POST_PGGB_DIR/merged_graph.gfa" -S -t "$THREADS" > "$POST_PGGB_DIR/merged_graphs.info.tsv"
echo "Done generatig graph statistics!"
echo "------------------------------------------------"

# ---------------------------------------------- #
echo "Getting sample list for ordered growth statistics..."
awk '$1=="group" && $3=="bp" {print $2}'  "$POST_PGGB_DIR/merged_graphs.info.tsv" >  "$POST_PGGB_DIR/sample_order.txt"
NUM_SAMPLES=$(wc -l < $POST_PGGB_DIR/sample_order.txt )
echo "Done! Number of samples: $NUM_SAMPLES in pangenome graph."
echo "------------------------------------------------"

# Panacus: Ordered histgrowth
echo "Generating ordered growth table..."
RUST_LOG=info panacus ordered-histgrowth -c bp -t "$THREADS" -l 1,2,$NUM_SAMPLES -S -O "$POST_PGGB_DIR/sample_order.txt" "$POST_PGGB_DIR/merged_graph.gfa" > "$POST_PGGB_DIR/ordered-histgrowth.bp.tsv"
RUST_LOG=info panacus ordered-histgrowth -c node -t "$THREADS" -l 1,2,$NUM_SAMPLES -S -O "$POST_PGGB_DIR/sample_order.txt" "$POST_PGGB_DIR/merged_graph.gfa" > "$POST_PGGB_DIR/ordered-histgrowth.node.tsv"
echo "Done generating ordered growth table!"
echo "------------------------------------------------"


# Panacus: Visualize graph growth
echo "Generating ordered growth visualization..."
panacus-visualize "$POST_PGGB_DIR/ordered-histgrowth.bp.tsv" > "$POST_PGGB_DIR/ordered-histgrowth.bp.pdf"
panacus-visualize "$POST_PGGB_DIR/ordered-histgrowth.node.tsv" > "$POST_PGGB_DIR/ordered-histgrowth.node.pdf"
echo "Done! Visualization in PDF written to: $POST_PGGB_DIR/ordered-histgrowth.bp.pdf"
echo "------------------------------------------------"

# Panacus group by node, presence-absence of nodes
echo "Generating panacus table..."
panacus table --count node --groupby-sample --threads $THREADS  "$POST_PGGB_DIR/merged_graph.gfa" > "$POST_PGGB_DIR/merged_graphs.table.tsv"
echo "Done panacus generating table!"
echo "------------------------------------------------"

#
echo "Calculating Sum column"
awk 'NR==1 {print $0 "\tSum"; next} {sum=0; for(i=2;i<=NF;i++) sum+=$i; print $0 "\t" sum}' "$POST_PGGB_DIR/merged_graphs.table.tsv" > "$POST_PGGB_DIR/merged_graphs.table.sum.tsv"
echo "Done calculating sum"
echo "------------------------------------------------"

# quantify node coverage total
echo "Quantifying total node coverage..."
seq 0 $NUM_SAMPLES | while read i; do 
awk -v i="$i" '$2 == i' "$POST_PGGB_DIR/merged_graphs.table.sum.tsv" | wc -l; 
done
echo "Done quantifying node coverage!"
echo "------------------------------------------------"

echo "Categorizing nodes into core, private, and dispensable"
awk '$NF == $NUM_SAMPLES' "$POST_PGGB_DIR/merged_graphs.table.sum.tsv" | cut -f1 > "$POST_PGGB_DIR/core.nodes.tsv"
CORE_NODES=$(wc -l < $POST_PGGB_DIR/core.nodes.tsv)
echo "Total core nodes: $CORE_NODES"
echo "------------------------------------------------"

awk '$NF == 1' "$POST_PGGB_DIR/merged_graphs.table.sum.tsv" | cut -f1 > "$POST_PGGB_DIR/private.nodes.tsv"
PRIVATE_NODES=$(wc -l < $POST_PGGB_DIR/private.nodes.tsv)
echo "Total private nodes: $PRIVATE_NODES"
echo "------------------------------------------------"

awk '$NF> 1 && $NF < $NUM_SAMPLES' "$POST_PGGB_DIR/merged_graphs.table.sum.tsv" | cut -f1 > "$POST_PGGB_DIR/dispensable.nodes.tsv"
DISPENSABLE_NODES=$(wc -l < $POST_PGGB_DIR/dispensable.nodes.tsv)
echo "Total dispensable nodes: $DISPENSABLE_NODES"
echo "Done categorizing nodes!"
echo "------------------------------------------------"

} > "$LOG_DIR/post_pggb_${TIMESTAMP}.log" 2>&1
