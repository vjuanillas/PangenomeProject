library(data.table)
library(ggplot2)
library(parallel)
library(dplyr)

# set wd
setwd("/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p90")

##################
## Notes #########
##################
# Check the dataframe for obtaining genotypes. Proper index for column is needed for genotypes 
# Check combinations
# Implement SNPs soon

#####################
### Load VCF File ###
#####################

vcf <- fread("pggb_merged_final_clean.vcf")

# split by chromosome
vcf_raw <- split(vcf, by = "#CHROM")
print("Raw VCF loaded!")

vcf_info_12chrs <- rbindlist(lapply(1:length(vcf_raw), function(x) {
  vcf_raw_now <- vcf_raw[[x]]
  chr_now <- gsub("O_mer#0#out_", "", names(vcf_raw[x]))
  
  ##############################
  ### Get Allele information ###
  ##############################
  
  
  N_allele <- vcf_raw_now[,.(all=paste(REF,ALT,sep = ","))][,lapply(.SD, function(x) sapply(strsplit(x,","),function(y) length(y)))][,setnames(.SD,"N_allele")]
  N_bp_ref <- vcf_raw_now[, .(N_bp_ref = nchar(REF))]
  N_bp_diff <- vcf_raw_now[, .(all = paste(REF, ALT, sep = ","))][, lapply(.SD, function(x) sapply(strsplit(x, ","), 
                                                                                                function(y) max(nchar(y) - min(nchar(y)))))][, setnames(.SD, "N_bp_diff")]
  vcf_info <- vcf_raw_now[,c(1,2,4,5)][,setnames(.SD,old="#CHROM",new="CHROM")][,CHROM:=gsub("O_mer#0#out_","",CHROM)]
  vcf_info <- data.table(vcf_info, N_allele, N_bp_ref, N_bp_diff)
  
  ######################################
  ### Determine General Variant Type ###
  ######################################
  # Classify variants based on allele length differences
  vcf_info[N_bp_diff == 0 & N_bp_ref == 1, type := "SNP"]
  vcf_info[N_bp_diff == 0 & N_bp_ref < 50 & N_bp_ref != 1, type := "MNP"]
  vcf_info[N_bp_diff == 0 & N_bp_ref >= 50, type := "SV_complex"]
  vcf_info[N_bp_diff > 0 & N_bp_diff < 50, type := "INDEL"]
  vcf_info[N_bp_diff > 0 & N_bp_diff >= 50, type := "SV"]
  
  # Check the counts of each variant type
  vcf_info[, .N, by = type]
  
  # Assign variant length (N_bp)
  vcf_info[type == "SNP", N_bp := 1]
  vcf_info[type != "SNP", N_bp := N_bp_diff]
  
  ######################
  ### Get Genotype data
  ######################
  
  # Get genotype
  geno <- vcf_raw_now[, 10:15]
  print(paste("Genotype for the following obtained:",colnames(geno)))
  
  # combine tables
  vcf_info <- data.table(vcf_info, geno)
  
  # calculate missing data
  missing <- vcf_raw_now[, 10:15][,rowSums(.SD == ".")]
  
  # add data
  vcf_info <- data.table(vcf_info, missing)
  print(paste("Processing missing genotypes for chromosome done",chr_now))
  
  
  ##################################
  ### Determine Variant Subtypes ###
  ##################################
  # Classify SNPs and MNPs
  vcf_info[type == "SNP", subtype := "SNP"]
  vcf_info[type == "MNP", subtype := "MNP"]
  
  # Classify INDELs (2-49 bp), 
  vcf_info[type == "INDEL" & nchar(REF) == 1, subtype := "INS"]
  vcf_info[type == "INDEL" & nchar(ALT) == 1, subtype := "DEL"]
  vcf_info[type == "INDEL" & is.na(subtype), subtype := "INDEL_complex"]
  vcf_info[type == "INDEL" & N_allele > 2, subtype := "Multiallelic_INDEL"]
  
  # Classify SVs (>= 50 bp),
  vcf_info[type == "SV" & nchar(REF) == 1, subtype := "SV_INS"]
  vcf_info[type == "SV" & nchar(ALT) == 1, subtype := "SV_DEL"]
  vcf_info[type == "SV" & is.na(subtype), subtype := "SV_complex"]
  vcf_info[type == "SV_complex", subtype := "SV_complex"]
  vcf_info[type == "SV" & N_allele > 2, subtype := "Multiallelic_SV"]
  
  vcf_info[N_allele == 2, bi_multi := "Biallelic"]
  vcf_info[N_allele > 2, bi_multi := "Multiallelic"]
  
  # Check counts by type and subtype
  vcf_info[, .N, by = type]
  vcf_info[, .N, by = subtype]
  
  # save per-chromosome summary
  fwrite(vcf_info[order(POS)], file=paste0("VCF_summary_",chr_now,".txt"), sep="\t")
  
  vcf_info
  }))

# Save combined VCF information for all chromosomes
fwrite(vcf_info_12chrs,"VCF_INFO_allchrs.txt",sep = "\t",quote = F,row.names = F,col.names = T)

# making sure data frame exists
vcf_info_12chrs <- vcf_info_12chrs

# summarize variant number and lengths
INS_all <- vcf_info_12chrs %>% group_by(CHROM) %>% 
  filter(subtype == "INS") %>% 
  summarise(N_bubbles = n(), Total_length = sum(nchar(ALT)), Avg_length = Total_length/N_bubbles,.groups = "drop")
DEL_all <- vcf_info_12chrs %>% group_by(CHROM) %>% 
  filter(subtype == "DEL") %>% 
  summarise(N_bubbles = n(), Total_length = sum(nchar(REF)), Avg_length = Total_length/N_bubbles,.groups = "drop")
SV_INS_all <- vcf_info_12chrs %>% group_by(CHROM) %>% 
  filter(subtype == "SV_INS") %>% 
  summarise(N_bubbles = n(), Total_length = sum(nchar(ALT)), Avg_length = Total_length/N_bubbles,.groups = "drop")
SV_DEL_all <- vcf_info_12chrs %>% group_by(CHROM) %>% 
  filter(subtype == "SV_DEL") %>% 
  summarise(N_bubbles = n(), Total_length = sum(nchar(REF)), Avg_length = Total_length/N_bubbles,.groups = "drop")

# Save merged summaries
fwrite(INS_all, "allchrs_INS_summary.txt", sep = "\t", quote = F, row.names = F, col.names = T)
fwrite(DEL_all, "allchrs_DEL_summary.txt", sep = "\t", quote = F, row.names = F, col.names = T)
fwrite(SV_INS_all, "allchrs_SV_INS_summary.txt", sep = "\t", quote = F, row.names = F, col.names = T)
fwrite(SV_DEL_all, "allchrs_SV_DEL_summary.txt", sep = "\t", quote = F, row.names = F, col.names = T)

print("All per-chromosome files and merged summaries saved!")

#####################################
##### Summary of Variant Counts #####
#####################################

# split by chromosome clean table
vcf_clean <- split(vcf_info_12chrs, by = "CHROM")
print("Clean VCF loaded for gain/loss analysis!")

vcf_gene_12chrs <- lapply(1:length(vcf_clean), function(x) {
  vcf_gene_now <- vcf_clean[[x]]
  chr_now <- names(vcf_raw[x])
  
  # handle sites with complete data only
  pattern_var <- vcf_gene_now %>% filter(missing == 0)
  
  # collapse accession names on new column for pattern matching
  pattern_var$REF_accn <- apply(pattern_var[,10:15], 1, function(x) paste(names(x)[which(x == 0)], collapse = ","))
  pattern_var$ALT_accn <- apply(pattern_var[,10:15], 1, function(x) paste(names(x)[which(x == 1)], collapse = ","))
  
  # SNPs 
  # implement soon
  
  # INDELS
  pattern_var_INDEL_gain <- pattern_var %>% filter(subtype == "INS") %>% group_by(ALT_accn) %>% summarize(N_bubbles=n(), alt_length = sum(nchar(ALT)))
  pattern_var_INDEL_loss <- pattern_var %>% filter(subtype == "DEL") %>% group_by(ALT_accn) %>% summarize(N_bubbles=n(), ref_length = sum(nchar(REF)))
  
  # SV
  pattern_var_SV_gain <- pattern_var %>% filter(subtype == "SV_INS") %>% group_by(ALT_accn) %>% summarize(N_bubbles=n(), alt_length = sum(nchar(ALT)))
  pattern_var_SV_loss <- pattern_var %>% filter(subtype == "SV_DEL") %>% group_by(ALT_accn) %>% summarize(N_bubbles=n(), ref_length = sum(nchar(REF)))
  
  #####################
  ### combinations ####
  #####################
  # this will vary case to case
  singletons <- c("Ob","Og","Osi","Or","Osj","Osj(IRGSP)")
  lineage <- c("Osj,Osj(IRGSP)","Ob,Og","Or,Osj,Osj(IRGSP)","Or,Osi,Osj,Osj(IRGSP)")
  roots <- c("Ob,Og,Or,Osi,Osj,Osj(IRGSP)")
  
  combinations <- c(singletons, lineage, roots)
  
  # function to analyze combinations
  match_combination <- function(df, combinations) df %>% filter(ALT_accn %in% combinations)
  # summarizing INDEL parsimonious combinations
  pattern_var_INDEL_summary_gain <- match_combination(pattern_var_INDEL_gain,combinations)
  pattern_var_INDEL_summary_loss <- match_combination(pattern_var_INDEL_loss,combinations)
  
  # summarizing SV parsimonious combinations
  pattern_var_SV_summary_gain <- match_combination(pattern_var_SV_gain,combinations)
  pattern_var_SV_summary_loss <- match_combination(pattern_var_SV_loss,combinations)

  # save per-chromosome summary
  fwrite(pattern_var_INDEL_summary_gain, file=paste0("GAIN_INDEL_summary_",chr_now,".txt"), sep="\t")
  fwrite(pattern_var_INDEL_summary_loss, file=paste0("LOSS_INDEL_summary_",chr_now,".txt"), sep="\t")
  fwrite(pattern_var_SV_summary_gain, file=paste0("GAIN_SV_summary_",chr_now,".txt"), sep="\t")
  fwrite(pattern_var_SV_summary_loss, file=paste0("LOSS_SV_summary_",chr_now,".txt"), sep="\t")
  # Return gain and loss separately for INDELs and SVs
  return(list(
    gain_INDEL = pattern_var_INDEL_summary_gain,
    loss_INDEL = pattern_var_INDEL_summary_loss,
    gain_SV = pattern_var_SV_summary_gain,
    loss_SV = pattern_var_SV_summary_loss
  ))
})

# Separate gains and losses
gains_INDEL_all <- rbindlist(lapply(vcf_gene_12chrs, `[[`, "gain_INDEL"), fill = TRUE)
losses_INDEL_all <- rbindlist(lapply(vcf_gene_12chrs, `[[`, "loss_INDEL"), fill = TRUE)
gains_SV_all <- rbindlist(lapply(vcf_gene_12chrs, `[[`, "gain_SV"), fill = TRUE)
losses_SV_all <- rbindlist(lapply(vcf_gene_12chrs, `[[`, "loss_SV"), fill = TRUE)

# summarize
gains_INDEL_all <- gains_INDEL_all %>% group_by(ALT_accn) %>% summarise(N_bubbles = sum(N_bubbles),
                                                            gain_length = sum(alt_length),
                                                            .groups = "drop")
losses_INDEL_all <- losses_INDEL_all %>% group_by(ALT_accn) %>% summarise(N_bubbles = sum(N_bubbles),
                                                                  loss_length = sum(ref_length),
                                                                  .groups = "drop")
gains_SV_all <- gains_SV_all %>% group_by(ALT_accn) %>% summarise(N_bubbles = sum(N_bubbles),
                                                                        gain_length = sum(alt_length),
                                                                        .groups = "drop")
losses_SV_all <- losses_SV_all %>% group_by(ALT_accn) %>% summarise(N_bubbles = sum(N_bubbles),
                                                                  loss_length = sum(ref_length),
                                                                  .groups = "drop")

# Save merged summaries
fwrite(gains_INDEL_all, "allchrs_GAIN_INDEL_INFO.txt", sep = "\t", quote = F, row.names = F, col.names = T)
fwrite(losses_INDEL_all, "allchrs_LOSS_INDEL_INFO.txt", sep = "\t", quote = F, row.names = F, col.names = T)
fwrite(gains_SV_all, "allchrs_GAIN_SV_INFO.txt", sep = "\t", quote = F, row.names = F, col.names = T)
fwrite(losses_SV_all, "allchrs_LOSS_SV_INFO.txt", sep = "\t", quote = F, row.names = F, col.names = T)

print("All per-chromosome files and merged summaries saved!")