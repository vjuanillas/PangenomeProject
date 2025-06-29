#!/usr/bin/env python

"""

Extract the biallelic SV from bubble in the graphs
and output the mutation types. Modified version of Crysnanto et al. script on SV annotation

"""
import argparse
from dataclasses import dataclass

@dataclass
class svid:
    """
    An object to store the sv information
    """
    chromo: str
    pos: int
    svtype: str
    reflen: int
    nonreflen: int
    sourcenode: str
    bubref: str
    bubnonref: str
    sinknode: str


def sv_det(line):
    """
    Determines the SV types per line
    Input: Each line of the gfatools bubble output
    Output: Insertions or Deletions SV based on comparison between ref and non-ref allele
    """
    # sample line: chr01 11093 11118 4 2 s1,s191853,s2,s3
    chromo, pos, end, nodes, paths, nodelist = line.strip().split() # added end
    bubble = nodelist.split(",")[1:-1] # take nodes except first(source) and last(sink)
    sourcenode = nodelist.split(",")[0]
    sinknode = nodelist.split(",")[-1]
    if nodes in ["3","4"]:
        if nodes == "3": #either Insertion or Deletion
            rrank, nodelen = nodeinf[bubble[0]] # nodeinf refer to above cell
            if rrank > 0: # if rank of the first node inside the bubble is not the reference rank (means reference has len 0)
                svtype = "Insertion"
                reflen = 0
                nonreflen = nodelen
                bubnonref = bubble[0]
                bubref = 0 # 0 because no reference node
            else: # rank is on reference (means non-ref path has len 0)
                svtype = "Deletion"
                reflen = nodelen
                nonreflen = 0
                bubnonref = 0
                bubref = bubble[0]
        elif nodes == "4":
            if bubble[0] == bubble[1]: # typical cases are inversions so do not loop around the bubble
                return None
            for bub in bubble:
                rrank, nodelen = nodeinf[bub]
                if rrank > 0:
                    nonreflen = nodelen
                    bubnonref = bub
                else:
                    reflen = nodelen
                    bubref = bub
            svtype = "AltDel" if nonreflen < reflen else "AltIns"
        return chromo, pos, end, svtype, reflen, nonreflen, sourcenode, bubref, bubnonref, sinknode # added end
    else:
        return None # for biallelic paths with several nodes ; in this analysis, they are usually inversions

if __name__ == "__main__":

    # parse the assembly
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--combined", help="Combined Coverage file")
    parser.add_argument("-b", "--biallelic", help = "Biallelic Bubble tsv file")
    args = parser.parse_args()
    combined_coverage = args.combined
    biallelic_bubble = args.biallelic

    # container of the node including len and rank
    nodeinf = {}


    # s2 1389 1 348029 0
    # note, we are using our own script outputs. Revised for combined coverage output

    with open(combined_coverage) as infile:
        next(infile)
        for line in infile:
            line_comp = line.strip().split()
            nodeid = line_comp[0]
            nodelen = line_comp[1]
            chromo = line_comp[2]
            pos = line_comp[3]
            rrank = line_comp[4]
            nodeinf[nodeid] = [int(rrank), int(nodelen)] # key is nodeid, values are rank and the length

    with open(biallelic_bubble) as infile:
        for line in infile:
            if sv_det(line):
                print(*sv_det(line)) # it'll be better if the end coordinate is also included # added 08.17.2024

