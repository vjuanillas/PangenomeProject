# Gain and Loss Analysis
The gain and loss analysis begins with processing and deconstructing the variants from the pangenome. 

Make sure your current environment as the following programs installed:
1. panacus v.0.2.5
2. odgi v0.9.0-0-g1895f496
3. seqkit v2.6.0
4. RepeatMasker v4.1.7
5. Diamond v2.1.10

## Steps:
1. Pre-process the VCF file. The `Scripts/vcf_handle_v2.sh` should deconstruct, decompose, and normalize the previously generated vcf files. To do this, place all address of the vcf decomposed files such as seen in `vcf_gz.txt`. Then, execute `Scripts/vcf_handle_v2.sh`
```
bash vcf_handle_v2.sh
```
This should generate a cleaned and merged vcf file from each community. 
2. Next, we classify and filter the variants into classes: INDELs and SVs. The definition of INDELs is variants with length less than 50 bp. Meanwhile, variants with length greater than 50 bp are defined as SVs. Use and edit `Scripts/1_Process_vcf.R` in Rstudio to tail to your needs. There are sections that you need to edit inside the R script. This section is entitled:
`Per-taxon parsimonious combination`

Edit this part based on your parsimonious combination. 

```
Rscript 1_Process_vcf.R
```

You will obtain the following outputs:
| Name | Description |
| ---- | ----------- |
| whole_genome/VCF_RAW_allchrs.txt | Combined VCF information for all chromosomes including genotypes |
| whole_genome/VCF_stats_raw_summary.txt | Statistics per chromosome |
| whole_genome/VCF_stats_WG_raw_summary.txt | Statistics for whole genome |
| whole_genome/VCF_clean_variant_stats_perChrom.txt | Raw variant counts per chromosome | 
| whole_genome/VCF_clean_variant_stats_WG.txt | Raw variant counts for whole genome |
| whole_genome/VCF_CLEAN_allchrs.txt | Clean counts for all chromosomes |
| whole_genome/VCF_Parsimony_Not_raw_statistics.txt | Includes parsimonious and non-parsimonious counts |
| whole_genome/VCF_Parsimony_stats_perChrom.txt | Per chromosome VCF parsimony statistics |
| whole_genome/VCF_Parsimony_stats_perWG.txt | Whole genome parsimonious statistics |
| whole_genome/Gain_parsimony_all_summary.txt | Gain parsimony summary |
| whole_genome/Loss_parsimony_all_summary.txt | Loss parsimony summary |
| whole_genome/VCF_Clean_parsimony_labelled_gain_loss_final.txt | Cleaned parsimony statistics |
| whole_genome/GAIN_statistics_final.txt | GAIN statistics final |
| whole_genome/LOSS_statistics_final.txt | Loss statistics final |
whole_genome/VCF_clean_gain_forMAP.txt | gain with accession |
| whole_genome/VCF_clean_loss_forMAP.txt | loss with accessions |
|

3. Visualize the gain and loss using `Scripts/2_tree_viz_gainloss_v4.R`. Adjust the arrangement and the combinations to your context. 
4. To quantify the repeats for each variation included in the parsimony tree, extract the sequences using a custom script called `Scripts/Extract_fasta_Annotate_GAIN_LOSS.R`. Then, use the script `Scripts/3_Repeats_extract_visualize_repeats.R` to visualize the fractions. 
5. Deconstruct using IRGSP
6. bedtools intersect
7. Use the `Scripts/4_VCF_intersection.R` to analyze the genes for 