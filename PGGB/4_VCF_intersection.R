library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

setwd("/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p90_v2/gene_intersection")
vcf_intersection <- fread("VCF_to_OmerGFF_noChromosome.tsv")

colnames(vcf_intersection) <- c("CHR","POS_1","POS_2","NODE","REF","ALT","TYPE","N_bp","SUBTYPE","REF_accn","ALT_accn",
                                "CHR_gff3","START","END","gene_struc","NAME")
vcf_intersection_indels <- vcf_intersection %>% filter(SUBTYPE == c("SV_INS","SV_DEL","INS","DEL"))
# if "." = intergenic
vcf_intersection_indels$gene_struc <- ifelse(vcf_intersection_indels$gene_struc == ".", "Intergenic", vcf_intersection_indels$gene_struc)


######################
### Prioritization ###
######################

priority <- c("CDS"=1, "three_prime_UTR"=2, "five_prime_UTR"=3, "exon"=4, "mRNA"=5, "gene"=6, "Intergenic"=7)

# Assign priority and fallback for unknowns
vcf_intersect_nr <- vcf_intersection_indels %>%
  mutate(priority_rank = ifelse(gene_struc %in% names(priority), priority[gene_struc], 99)) %>%
  group_by(CHR, POS_1) %>%
  slice_min(order_by = priority_rank, n = 1, with_ties = FALSE) %>%
  ungroup()

vcf_intersect_nr$TRANSCRIPT <- str_extract(vcf_intersect_nr$NAME, "(?<=transcript:)[^;]+")

# summarize raw variant counts
vcf_intersect_nr <- data.table(vcf_intersect_nr) %>% filter(priority_rank != 99)
vcf_intersect_raw_stats <- vcf_intersect_nr[, .(
  N_bubble = .N,
  Variant_length = sum(N_bp),
  Average_length = sum(N_bp)/.N
), by = .(TYPE, SUBTYPE, gene_struc)]

# for export purposes
# vcf_intersect_nr_ <- vcf_intersect_nr %>% select(-priority_rank, -NAME) 
# fwrite(vcf_clean_raw_stats, "whole_genome/VCF_raw_variant_stats_perChrom.txt", sep = "\t", quote=F, row.names = F, col.names=T)

vcf_intersect_raw_stats_INDEL <- vcf_intersect_raw_stats %>% filter(TYPE == "INDEL")
vcf_intersect_raw_stats_SV <- vcf_intersect_raw_stats %>% filter(TYPE == "SV")
#####################
### Visualization ###
#####################
library(ggplot2)

# number of occurence
plot_gene_struc <- function(df, df_y, x_label="Subtype", y_label){
  ggplot(df, aes(x = SUBTYPE, y = {{df_y}}, fill = gene_struc)) +
    geom_bar(stat = "identity", colour = "black", size=0.2) +
    labs(
      x = "Subtype",
      y = y_label,
      fill = "Gene Structure",
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_classic() +
    theme(
      text = element_text(size = 12)
    )
}
a_indel <- plot_gene_struc(vcf_intersect_raw_stats_INDEL,N_bubble, x_label="Subtype", y="Number of Occurences")
b_indel <- plot_gene_struc(vcf_intersect_raw_stats_SV,N_bubble, x_label="Subtype", y="Number of Occurences")

# variant length
plot_Length <- function(df, df_y, x_label="Gene Structure", y_label="Value", y_scale=1){
  ggplot(df, aes(x = gene_struc, y = {{df_y}}/y_scale, fill = SUBTYPE)) +
    geom_bar(stat = "identity", position = "dodge", colour = "black", linewidth = 0.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = x_label, 
         y = y_label, 
         fill = "Type") +
    theme_classic() +
    theme(strip.background = element_blank(), 
          strip.text = element_text(size = 12, face = "bold"))
}
plot_Length(vcf_intersect_raw_stats_INDEL, Variant_length, x_label="Gene Structure","Cumulative length (in Kp)",1e3)
plot_Length(vcf_intersect_raw_stats_SV, Variant_length, x_label="Gene Structure","Cumulative length (in Mb)", 1e6)

# average length
plot_Length(vcf_intersect_raw_stats_INDEL, Average_length, x_label="Gene Structure","Average length (in bp)")
plot_Length(vcf_intersect_raw_stats_SV, Average_length, x_label="Gene Structure","Average length (in bp)")

## convert mRNA and gene to intron??
## appropriate scaling in average lengths