#!/usr/bin/env python

"""
Parse all GFA files with dc tag and output a data frame with combined coverage per segment and edge
"""

import argparse
import pandas as pd 
import os 

def parse_args():
	parser = argparse.ArgumentParser(description = "Combine GFA files with dc tags")
	parser.add_argument("-g", "--graph", help = "Pangenome Graph name")
	parser.add_argument("-a", "--assembly", help = "Assembly Names", nargs="+")
	parser.add_argument("-d", "--directory", help = "Directory of GFA files with dc tags")
	parser.add_argument("-o", "--output", help = "Directory for output files")
	return parser.parse_args()

"""
In each gfa file with dc tag, the following format is observed:
S/L segment_ID sequence segment_length chromosome offset rank coverage
S 	s3 	CGTG 	LN:i:1 	SN:Z:1 	SO:i:1202 	SR:i:0 	dc:f:0

This function should parse the gfa files
Input: GFA files with dc tag per assembly
Output: a tuple of segment_ID, segment_length, start_chromosome, offset, rank, and coverage
"""
def parser_segment_coverage(line):
	line_comp = line.strip().split() # strip all whitespaces and split into components
	segment_id = line_comp[1] # take segment ID on the first tag
	segment_length = len(line_comp[2]) # take length of the sequence
	start_chromosome = line_comp[4].split(":")[2] # take chromosome tag, split by colon, then take the third element
	offset = line_comp[5].split(":")[2] # take offset, split at colon, and take the third element
	segment_rank = line_comp[-2].split(":")[2] # take rank tags, split by colon, and take the third element
	coverage = line_comp[-1].split(":")[2] # take dc tags, split at colon, and take the third element

	return segment_id, segment_length, start_chromosome, offset, segment_rank, coverage

def parser_edge_coverage(line):
	line_comp = line.strip().split() 
	if line_comp[2] == "-" or line_comp[4] == "-": # reverse versus forward direction
		parent = line_comp[3]
		child = line_comp[1]
	else:
		parent = line_comp[1]
		child = line_comp[3]
	coverage = int(line_comp[-1].split(":")[2]) # at the last part of line, the dc tag, split and take the integer at [2]
	coverage_rep = 1 if coverage > 0 else 0
	return (parent, child, coverage_rep)

if __name__ == "__main__":
	args = parse_args() # parse the arguments
	graph_name = args.graph # assign variables to each arguments
	asms = args.assembly # names of the assemblies
	folder = args.directory # name of the input directory
	out = args.output # name of the output directory

	for asm in asms: # loop through each assembly names 
		with open(f"{folder}/{asm}_{graph_name}.gfa") as infile:
			if asm == asms[0]:
			# for the first assembly, parse all columns
			# initialize data frame
				combined = pd.DataFrame([parser_segment_coverage(line) for line in infile if line.startswith("S")],
					columns=["SegmentId", "SegmentLen", "StartChrom", "Offset", "SegRank", asm])
			else:
				# for subsequent assemblies, only take the coverage at the last column
				append_cov = pd.DataFrame([[parser_segment_coverage(line)[0], parser_segment_coverage(line)[-1]] for line in infile if line.startswith("S")],
					columns =["SegmentId", asm])
				# merge the existing data frame with the new data frame based on the segment ID and add on the outer column
				combined = pd.merge(combined, append_cov, on="SegmentId", how="outer")

		with open(f"{folder}/{asm}_{graph_name}.gfa") as infile:
			if asm == asms[0]:
				combined_edge = pd.DataFrame([parser_edge_coverage(line) for line in infile if line.startswith("L")],
					columns=["parent_segment","child_segment", asm])
			else:
				append_edge = pd.DataFrame([parser_edge_coverage(line) for line in infile if line.startswith("L")],
					columns=["parent_segment", "child_segment", asm])
				combined_edge = pd.merge(combined_edge, append_edge, on=["parent_segment", "child_segment"], how="outer")




combined.fillna(0) # fill all NAs with 0
# save file as csv
if not os.path.exists(out): # if the directory does not exist, create the directory
	os.makedirs(out)

combined.to_csv(f"{out}/01_combined_coverage.tsv", sep=" ", index=False)
combined_edge.to_csv(f"{out}/01_combined_coverage_edge.tsv", sep=" ", index=False)

		





