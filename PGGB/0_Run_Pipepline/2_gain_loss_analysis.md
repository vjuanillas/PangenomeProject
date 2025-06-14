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

2. 