# Processing Minigraph's output 

1. 01_coverage.py
Description:
Parses the GFA files with the dc tag and combine for all assemblies specified. 
This script follows the naming format of a file {asm}_{graph}.gfa. When executing the command, the name of the input directory is to be specified. This directory should contain all the gfa files to be parsed, following the naming format. Additionally, an output directory name is to be passed where the files will be saved. If the file is not yet existing, the python script will create the directory. 

INPUT: Graph name, assembly names, input directory, and output directory
OUTPUT: 01_combined_coverage.tsv, a file with all the coverages parsed from the gfa remapping
```
python 01_coverage.py -g nip -a IRGSP nh232 cw02 nh236 nh286 nh273 -d remapping_files -o coverage_v2
```
# update 07.27 - added parser edge coverage as input in structural variations calling

2. 02_Presence_absence.R
Description:
Creates a presence-absence matrix for the segments based on the combined coverages. There are two conditions to be fullfilled for a segment to be present when the assemblies are remapped back to the graph (x = 1). The conditions are the following: (1) The coverage of the assembly when mapped back to the graph on the particular segment is greater than 0 (dc > 0) and (2) the segment's rank and the assigned rank of the assembly by minigraph is the same (assembly rank == segment rank). Note that the script is to be called on the directory where the file will be generated unlike the previous script. 

INPUT: the combined_coverage file from previous script; 01_combined_coverage.tsv
OUTPUT:  
1. 02_Presence_Absence_ALL.tsv -  all of the information after parsing
2. 02_Presence_Absence_Matrix_only.tsv - the presence-absence matrix denoted by 1 and 0 but names are not included
3. 02_Presence_Absence_Names_only.tsv - the presence-absence matrix based on the accession names but not the values of 1 and 0
```
Rscript ../../Scripts/02_Segment_matrices.R 01_combined_coverage.tsv
```

3. 03_Ref_Analysis.R
Description:
Outputs several tsv values that contains the classification of the pangenome segments: core, dispensable, and private segments. This script contains three functions: 1) reference_categorization is a function that takes the previously constructed dataframe with presence-absence matrices from script 02 and outputs the categorization statistics such as the number of segments and the cumulative length of the categories; 2) intersection is a function that takes the previously constructed dataframe and outputs a list with the combinations necessary for the intersection plot; and lastly, 3) reference_sharing is a function that takes a similar input and quantifies the cumulative lengths of the different combinations of accessions and their corresponding segments. 

INPUT: matrix from script 02
OUTPUT: 
1. 03_Ref_Category_CountLen.tsv - summarizes the number of segments and their corresponding length per category (core, dispensable, and private)
2. 03_Ref_Category_Core.tsv, 03_Ref_Category_Dispensable.tsv, and 03_Ref_Category_Private.tsv - tsv files containing segment information including ranks and segment offsets for each category. 
3. segment_list - a list outputted after running intersection function which can be used as a direct input to create the intersection plot
4. plt - a plot object after running reference_sharing function which can be run in R to visualize the different cumulative lengths of segments per combination of accession in the pangenome based on the presence-absence matrix.
Example command: this script is to be sourced in R and its function to be used separately

4. 03.2_NonRef_Analysis.R
Description:
This script requires the matrix outputted by script #2. It filters the ranks > 0 to determine the nonreference sequences and outputs an image of the sharing patterns between the assemblies, including the reference genome. It needs the name of the reference genome for the determination of nonreference segments as it needs to filter columns in R.

INPUT: matrix from script 02
OUTPUT: 
1. 02_nonref_shared_count.png - returns an image of the sharing counts of nonreference sequences
2. 02_nonref_shared_len.png -  returns an image of the cumulative lengths of the shared nonreference sequences
3. 02_02_nonref_sharing_intersection - returns an intersection
```
Rscript ../../Scripts/03.2_NonRef_analysis.R 02_Presence_Absence_ALL.tsv "IRGSP"
```

5. 04_GO_viz.R
Description: Takes a dataframe and creates a combined plot for visualization of gene ontology terms for each specificed category and number of genes. It has two functions: 1) go_visualize_1 function is strictly for the three categories of pangenomes ; 2) go_visualize_2 is for any dataframe without strict segment labels
OUTPUT:
1. ggplot object - can be visualized using R print
Example command: this script is to be sourced in R and its function to be used separately

6. 05.1_get_bialsv.py (Note: with variation)
Description: This is a modified version of Crysnanto et al. script on bialsv annotation. Note that this is to be used right after the linux terminal commands are executed for the separation of bialleles and multialleles. It takes the combined coverage and biallelic tsv file and annotate the types of structural variations: 1) Insertion, 2) Deletion, 3) AltIns, and 4) AltDel. 
Input: Combined_coverage.tsv and asm5.nip.biallelic.bubble.tsv (from linux commands)
Output: asm5.nip.biallelic_sv.tsv - contains the biallelic sv columns plus the structural variation types
```
python 05_get_bialsv.py -b biallelic.bubble.tsv -c combined_coverage.tsv > asm5.nip.biallelic_sv.tsv
```

7. 05.2_trace_path.py
Description: Use the biallelic SV output from script 0x_get_bialsv.py. This script processes each of the bubble files generated by
minigraph --call and traces the path of each assemble through the bubbles. It reads a directory with a strict file name syntax {asm}.alleles.bed and {asm} becomes the name of this file inside a temporary dictionary. 
Input: output from 05.1 script, path to the bubble files, and output path
Output: test_05_biallelic_asm_paths.tsv - tsv file containing biallelic bubbles, its type, and the ref and alt paths and asms
```
python 0x_trace_path.py -b /Users/jongpaduhilao/Desktop/LAB_Files/Initial_Pangenome_analysis/Trial_4/crysnanto_bubble/asm5.nip.biallelic_sv.tsv -a /Users/jongpaduhilao/Desktop/LAB_Files/Initial_Pangenome_analysis/bubble/04_bubble -o test_05_biallelic_asm_paths
```
update 11.27 - hashed the chromosome handling due to redundancy and leading to indexing error

Use bedtools for intersection with annotation files. 

8. 06.1_annot_sv.py
Description: Uses the bedtools intersect output between traced_paths and the transcripts.gff of IRGSP annotation. It parses the intersection output and returns what type of features does the structural variations fall in. Particularly, the script determines whether the bubbles are in mRNA, exon, 5' UTR, 3' UTR, intergenic, and CDS.
Input: 05x4_biallelic_asm2annot.intersect.tsv
Output: 06_biallelic_sv_annot.tsv
```
python Scripts/06.1_annot_sv.py -i Raw_data/data/05x4_biallelic_asm2annot.intersect.tsv -o 06_biallelic_sv_annot.tsv
``` 
 

