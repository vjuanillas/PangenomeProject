library(data.table)
library(dplyr)
library(tidyverse)

# set wd
# change the path to vcf if needed

setwd("/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p90_v2")

#####################
### Load VCF File ###
#####################

vcf <- fread("pggb_merged_final_clean_v2.vcf")

# split by chromosome
vcf_raw <- split(vcf, by = "#CHROM")
print("Raw VCF loaded!")

# create directory for output inside work_dir

if (!dir.exists("per_chromosome")) {
  dir.create("per_chromosome")
}

if (!dir.exists("whole_genome")) {
  dir.create("whole_genome")
}

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
  vcf_info <- vcf_raw_now[,c(1,2,3,4,5)][,setnames(.SD,old="#CHROM",new="CHROM")][,CHROM:=gsub("O_mer#0#out_","",CHROM)]
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
  
  #########################
  ### Get Genotype data ###
  #########################
  
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

  # save per-chromosome processed VCF files
  fwrite(vcf_info[order(POS)], file=paste0("per_chromosome/VCF_RAW_",chr_now,".txt"), sep="\t")
  
  vcf_info
  }))

# Save combined VCF information for all chromosomes including genotypes
fwrite(vcf_info_12chrs,"whole_genome/VCF_RAW_allchrs.txt",sep = "\t",quote = F,row.names = F,col.names = T)

# making sure data frame exists
vcf_info_12chrs <- vcf_info_12chrs

# summarize general statistics for variants per type per chromosome
vcf_summary <- vcf_info_12chrs[, .(
  N_bubble = .N,
  MissingRows = sum(rowSums(.SD == ".") > 0),  # Count rows with at least one "."
  REF_length = sum(nchar(REF)),  # Total length of REF column
  ALT_length = sum(nchar(ALT))   # Total length of ALT column
), by = .(type, subtype, CHROM), .SDcols = 11:16] # adjust column indices

# save statistics
fwrite(vcf_summary, "whole_genome/VCF_stats_summary.txt", sep = "\t", quote = F, row.names = F, col.names = T)

print("All per-chromosome files and merged summaries saved!")

################################################
##### Summary of Variant Counts: Parsimony #####
################################################

# handle sites with complete data only
vcf_clean <- vcf_info_12chrs %>% filter(missing == 0)

# collapse accession names on new column for pattern matching
vcf_clean$REF_accn <- apply(vcf_clean[,11:16], 1, function(x) paste(names(x)[which(x == 0)], collapse = ","))
vcf_clean$ALT_accn <- apply(vcf_clean[,11:16], 1, function(x) paste(names(x)[which(x == 1)], collapse = ","))

# summarize raw variant counts
vcf_clean_raw_stats <- vcf_clean[, .(
  N_bubble = .N,
  REF_length = sum(nchar(REF)),  # Total length of REF column
  ALT_length = sum(nchar(ALT))   # Total length of ALT column
), by = .(CHROM, subtype)]

fwrite(vcf_clean_raw_stats, "whole_genome/VCF_raw_variant_stats_perChrom.txt", sep = "\t", quote=F, row.names = F, col.names=T)

vcf_clean_raw_stats <- vcf_clean[, .(
  N_bubble = .N,
  REF_length = sum(nchar(REF)),  # Total length of REF column
  ALT_length = sum(nchar(ALT))   # Total length of ALT column
), by = .(subtype)]

fwrite(vcf_clean_raw_stats, "whole_genome/VCF_raw_variant_stats_WG.txt", sep = "\t", quote=F, row.names = F, col.names=T)


#####################
### combinations ####
#####################
# this will vary case to case
singletons <- c("Ob","Og","Osi","Or","Osj","Osj(IRGSP)")
lineage <- c("Osj,Osj(IRGSP)","Ob,Og","Or,Osj,Osj(IRGSP)","Or,Osi,Osj,Osj(IRGSP)")
roots <- c("Ob,Og,Or,Osi,Osj,Osj(IRGSP)") # filter this out which is 1s in all for alternative

# update combinations
inverse_pattern <- c("Ob,Og,Or,Osi,Osj","Ob,Og,Or,Osi,Osj(IRGSP)","Ob,Og,Osi,Osj,Osj(IRGSP)",
                 "Ob,Og,Or,Osj,Osj(IRGSP)","Ob,Or,Osi,Osj,Osj(IRGSP)","Og,Or,Osi,Osj,Osj(IRGSP)",
                 "Ob,Og,Or,Osi","Ob,Og,Osi")

combinations <- c(singletons, lineage, inverse_pattern)

# label parsimonious combinations
vcf_clean$par_not <- ifelse(vcf_clean$ALT_accn %in% combinations, "par", "not")

# Save combined VCF information for all chromosomes with parsimony combinations
fwrite(vcf_clean,"whole_genome/VCF_CLEAN_allchrs.txt",sep = "\t",quote = F,row.names = F,col.names = T)

vcf_clean <- vcf_clean %>% filter(!ALT_accn %in% roots & ALT_accn != "") # filter out the root prior to quantification of parsimonious combinations

# summarize parsimonious from non-parsimonious
vcf_parsimony_stats <- vcf_clean[, .(
  N_bubble = .N,
  REF_length = sum(nchar(REF)),  # Total length of REF column
  ALT_length = sum(nchar(ALT))   # Total length of ALT column
), by = .(CHROM, par_not, subtype)]

fwrite(vcf_parsimony_stats, "whole_genome/VCF_Parsimony_stats.txt", sep = "\t", quote=F, row.names = F, col.names=T)

# summarize whole-genome tree value for gain and loss 
gain_InsSV <- vcf_clean[subtype %in% c("INS","SV_INS"),
                         .(N_bubble = .N,
                           GAIN_length = sum(N_bp)
                           ),
                         by = .(CHROM, subtype, ALT_accn, par_not)
                         ]
setorder(gain_InsSV, subtype)

loss_DelSV <- vcf_clean[subtype %in% c("DEL","SV_DEL"),
                         .(N_bubble = .N,
                           LOSS_length = sum(N_bp)
                         ),
                         by = .(CHROM, subtype, ALT_accn, par_not)
]

setorder(loss_DelSV, subtype)


fwrite(gain_InsSV, "whole_genome/Gain_parsimony_all_summary.txt", sep = "\t", quote=F, row.names = F, col.names=T)
fwrite(loss_DelSV, "whole_genome/Loss_parsimony_all_summary.txt", sep = "\t", quote=F, row.names = F, col.names=T)

#############################################
### Per-taxon parsimonious combination ###
#############################################

# Define taxon groups and corresponding insertion and deletion rules
# gain
gain_taxon_groups <- list(
  "Osj(IRGSP)" = list(Ins = "Osj(IRGSP)", Del = "Ob,Og,Or,Osi,Osj"),
  "Osj" = list(Ins = "Osj", Del = "Ob,Og,Or,Osi,Osj(IRGSP)"),
  "Or" = list(Ins = "Or", Del = "Ob,Og,Osi,Osj,Osj(IRGSP)"),
  "Osi" = list(Ins = "Osi", Del = "Ob,Og,Or,Osj,Osj(IRGSP)"),
  "Og" = list(Ins = "Og", Del = "Ob,Or,Osi,Osj,Osj(IRGSP)"),
  "Ob" = list(Ins = "Ob", Del = "Og,Or,Osi,Osj,Osj(IRGSP)"),
  "Osj,Osj(IRGSP)" = list(Ins = "Osj,Osj(IRGSP)", Del = "Ob,Og,Or,Osi"),
  "Or,Osj,Osj(IRGSP)" = list(Ins = "Or,Osj,Osj(IRGSP)", Del = "Ob,Og,Osi"),
  "Or,Osi,Osj,Osj(IRGSP)" = list(Ins = "Or,Osi,Osj,Osj(IRGSP)"),
  "Ob,Og" = list(Ins = "Ob,Og")
)

# loss
loss_taxon_groups <- list(
  "Osj(IRGSP)" = list(Del = "Osj(IRGSP)", Ins = "Ob,Og,Or,Osi,Osj"),
  "Osj" = list(Del = "Osj", Ins = "Ob,Og,Or,Osi,Osj(IRGSP)"),
  "Or" = list(Del = "Or", Ins = "Ob,Og,Osi,Osj,Osj(IRGSP)"),
  "Osi" = list(Del = "Osi", Ins = "Ob,Og,Or,Osj,Osj(IRGSP)"),
  "Og" = list(Del = "Og", Ins = "Ob,Or,Osi,Osj,Osj(IRGSP)"),
  "Ob" = list(Del = "Ob", Ins = "Og,Or,Osi,Osj,Osj(IRGSP)"),
  "Osj,Osj(IRGSP)" = list(Del = "Osj,Osj(IRGSP)", Ins = "Ob,Og,Or,Osi"),
  "Or,Osj,Osj(IRGSP)" = list(Del = "Or,Osj,Osj(IRGSP)", Ins = "Ob,Og,Osi"),
  "Or,Osi,Osj,Osj(IRGSP)" = list(Del = "Or,Osi,Osj,Osj(IRGSP)"),
  "Ob,Og" = list(Del = "Ob,Og")
)

# for annotation purposes
vcf_clean_labelled <- vcf_clean %>% filter(subtype %in% c("INS","DEL","SV_INS","SV_DEL")) %>% filter(par_not == "par") %>% mutate(Gain_Loss = NA, Taxon = NA)
# Loop over each taxon group and label gain or loss
for (group in names(gain_taxon_groups)) {
  ins_pattern <- gain_taxon_groups[[group]]$Ins
  del_pattern <- gain_taxon_groups[[group]]$Del
    
  vcf_clean_labelled <- vcf_clean_labelled %>% mutate(Gain_Loss = ifelse(subtype %in% c("INS","SV_INS") & ALT_accn %in% ins_pattern,"Gain",
                                                                         ifelse(subtype %in% c("DEL","SV_DEL") & ALT_accn %in% del_pattern,"Gain",Gain_Loss)),
                                                        Taxon = ifelse((subtype %in% c("INS", "SV_INS") & ALT_accn %in% ins_pattern) | 
                                                                         (subtype %in% c("DEL", "SV_DEL") & ALT_accn %in% del_pattern),
                                                                       group, Taxon)
  )
}

for (group in names(loss_taxon_groups)) {
  ins_pattern <- loss_taxon_groups[[group]]$Ins
  del_pattern <- loss_taxon_groups[[group]]$Del
  
  vcf_clean_labelled <- vcf_clean_labelled %>% mutate(Gain_Loss = ifelse(subtype %in% c("DEL","SV_DEL") & ALT_accn %in% del_pattern, "Loss",
                                                                         ifelse(subtype %in% c("INS","SV_INS") & ALT_accn %in% ins_pattern, "Loss", Gain_Loss)),
                                                      Taxon = ifelse((subtype %in% c("DEL", "SV_DEL") & ALT_accn %in% del_pattern) | 
                                                                       (subtype %in% c("INS", "SV_INS") & ALT_accn %in% ins_pattern),
                                                                     group, Taxon)
  )
}

## fwrite vcf_clean_labelled with ID for alt_accn 
# gain
vcf_clean_gain <- vcf_clean_labelled %>% filter(Gain_Loss == "Gain") 
vcf_clean_gain <- vcf_clean_gain[,.(N_bubble = .N, GAIN_length = sum(N_bp)), # ifelse DEL, take Ref as length of gain
                                 by = .(CHROM, POS, ID, REF, ALT, type, subtype, ALT_accn, Taxon)] %>% rename(Type = type)
setorder(vcf_clean_gain, subtype)

vcf_clean_gain_stats <- vcf_clean_gain %>% group_by(CHROM, Type, Taxon) %>% 
  summarize(Total_Gain_Nbubble = sum(N_bubble, na.rm = TRUE), Total_Gain_Length = sum(GAIN_length, na.rm = TRUE))

# loss
vcf_clean_loss <- vcf_clean_labelled %>% filter(Gain_Loss == "Loss")
vcf_clean_loss <- vcf_clean_loss[,.(N_bubble = .N, LOSS_length = sum(N_bp)),
                                 by = .(CHROM, POS, ID, REF, ALT, type, subtype, ALT_accn, Taxon)] %>% rename(Type = type)
setorder(vcf_clean_loss, subtype)

vcf_clean_loss_stats <- vcf_clean_loss %>% group_by(CHROM, Type, Taxon) %>%
  summarize(Total_Loss_Nbubble = sum(N_bubble, na.rm = TRUE), Total_Loss_Length = sum(LOSS_length, na.rm=TRUE))

# save VCF files
fwrite(vcf_clean_gain_stats, "whole_genome/GAIN_statistics_final.txt", sep = "\t", quote = F, row.names = F, col.names = T)
fwrite(vcf_clean_loss_stats, "whole_genome/LOSS_statistics_final.txt", sep = "\t", quote = F, row.names = F, col.names = T)
# for mapping coordinates
fwrite(vcf_clean_gain, "whole_genome/VCF_clean_gain_forMAP.txt", sep = "\t", quote = F, row.names = F, col.names = T)
fwrite(vcf_clean_loss, "whole_genome/VCF_clean_loss_forMAP.txt", sep = "\t", quote = F, row.names = F, col.names = T)

