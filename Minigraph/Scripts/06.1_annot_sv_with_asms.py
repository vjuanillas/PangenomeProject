#!/usr/bin/env python3

import csv
import sys
import argparse

def get_feature_table(mergedbreak):
    """
    Parse the TSV line to get the ID, feature type, SV type, reference_path, and alternative_path.
    Yields: [svid, feature_type, ID, sv_type, reference_path, alternative_path]
    """
    for feature in mergedbreak:
        try:
            svid = feature[3]  # SVID is in the 4th column
            sv_type = feature[4]  # type of SV
            reference_path = feature[7]  # reference_path is in the 6th column
            alternative_path = feature[8]  # alternative_path is in the 7th column
            feature_type = feature[11]  # mRNA and gene structures
            attributes = feature[-1].strip()
            
            # Parse attributes
            attr_dict = dict(x.split('=') for x in attributes.split(';') if '=' in x)
            
            # Get ID, defaulting to the 'Parent' if 'ID' is not present
            ID = attr_dict.get('ID', attr_dict.get('Parent', 'unknown'))
            # Convert transcript ID to locus ID
            locus_ID = ID.replace('t', 'g').split('-')[0] if 't' in ID else ID
            
            yield [svid, feature_type, locus_ID, sv_type, reference_path, alternative_path]
        except Exception as e:
            print(f"Error processing line: {feature}", file=sys.stderr)
            print(f"Error message: {str(e)}", file=sys.stderr)
            continue

def extract_important_feature(prevfeat, curfeat):
    """
    Return the most important feature based on priority.
    """
    priority = {"CDS": 6, "five_prime_UTR": 5, "three_prime_UTR": 4, "exon": 3, "mRNA": 2, "gene": 1}
    prevprior = priority.get(prevfeat, 0)
    curprior = priority.get(curfeat, 0)
    if curprior == 0 and prevprior == 0:
        return "intergenic"
    return curfeat if curprior >= prevprior else prevfeat

def process_features(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        tsv_reader = csv.reader(infile, delimiter='\t')
        feature_gen = get_feature_table(tsv_reader)
        
        current_svid = None
        current_type = "intergenic"
        current_id = ""
        current_sv_type = ""
        current_ref_path = ""
        current_alt_path = ""
        
        for svid, feat_type, feat_id, sv_type, ref_path, alt_path in feature_gen:
            if svid != current_svid:
                if current_svid:  # Write the previous entry
                    outfile.write(f"{current_svid}\t{current_sv_type}\t{current_ref_path}\t{current_alt_path}\t{current_type}\t{current_id}\n")
                current_svid = svid
                current_type = "intergenic"
                current_id = feat_id
                current_sv_type = sv_type
                current_ref_path = ref_path
                current_alt_path = alt_path
            
            current_type = extract_important_feature(current_type, feat_type)
            if current_type != "intergenic":
                current_id = feat_id
        
        # Write the last entry
        if current_svid:
            outfile.write(f"{current_svid}\t{current_sv_type}\t{current_ref_path}\t{current_alt_path}\t{current_type}\t{current_id}\n")

def main():
    parser = argparse.ArgumentParser(description="Process merged break file and annotate structural variations.")
    parser.add_argument('-i', '--input', required=True, help='Input TSV file')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    
    args = parser.parse_args()

    process_features(args.input, args.output)
    print(f"Annotation complete. Results written to {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()