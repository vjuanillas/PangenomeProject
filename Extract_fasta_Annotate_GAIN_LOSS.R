library(data.table)
library(dplyr)
library(stringr)

setwd("/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p90_v2")

vcf_merge_gain <- fread("whole_genome/VCF_clean_gain_forMAP.txt")
vcf_merge_loss <- fread("whole_genome/VCF_clean_loss_forMAP.txt")

# create dir

if (!dir.exists("fasta_annot")) {
  dir.create("fasta_annot")
}

if (!dir.exists("fasta_annot/indel")) {
  dir.create("fasta_annot/sv")
}

if (!dir.exists("fasta_annot/sv")) {
  dir.create("fasta_annot/sv")
}

# Step 1: Select relevant columns from vcf_gain
vcf_tsv_gain_seq <- vcf_merge_gain %>%
  select(CHROM, POS, ID, REF, ALT, Type, subtype, ALT_accn, Taxon)
vcf_tsv_loss_seq <- vcf_merge_loss %>%
  select(CHROM, POS, ID, REF, ALT, Type, subtype, ALT_accn, Taxon)

# Step 2: Create the HEADER column by collapsing CHROM, POS, and ID
vcf_tsv_gain_seq <- vcf_tsv_gain_seq %>%
  mutate(ID = gsub(">", "#", ID), HEADER = paste(CHROM, POS, ID, sep = "_"))
vcf_tsv_loss_seq <- vcf_tsv_loss_seq %>%
  mutate(ID = gsub(">", "#", ID), HEADER = paste(CHROM, POS, ID, sep = "_"))


# Step 3 & 4: Extract SEQ based on subtype and remove the first character
vcf_tsv_gain_seq <- vcf_tsv_gain_seq %>%
  mutate(SEQ = ifelse(subtype %in% c("INS", "SV_INS"), ALT, REF)) %>%
  mutate(SEQ = substr(SEQ, 2, nchar(SEQ))) %>%  # Remove first character
  select(-CHROM,-POS, -ID) %>% select(HEADER, SEQ, everything())
vcf_tsv_loss_seq <- vcf_tsv_loss_seq %>%
  mutate(SEQ = ifelse(subtype %in% c("INS", "SV_INS"), ALT, REF)) %>%
  mutate(SEQ = substr(SEQ, 2, nchar(SEQ))) %>%  # Remove first character
  select(-CHROM,-POS, -ID) %>% select(HEADER, SEQ, everything())

# split fasta files to respective taxons
taxon_split <- function(df,type,taxon) {
  taxon_split_name <- df %>% filter(Type == type & Taxon == taxon)
  return(taxon_split_name)
}

# Step 5 & 6: Write to a formatted txt file
write_fasta <- function(df, file_path) {
  con <- file(file_path, "w")  # Open connection for writing
  on.exit(close(con))  # Ensures connection closes even if an error occurs
  
  for (i in seq_len(nrow(df))) {
    header <- paste0(">", df$HEADER[i])  # Format header
    seq <- df$SEQ[i]
    
    # Break sequence into lines of max 60 characters
    seq_lines <- str_wrap(seq, width = 60, exdent = 0)
    
    # Write to file
    writeLines(c(header, seq_lines), con)
  }
}

taxon_split2fasta <- function(df, type, taxon, file_path) {
  taxon_split_now <- taxon_split(df, type, taxon)
  write_fasta(taxon_split_now, file_path)
}

######################
### IMPLEMENTATION ###
######################

###################
### GAIN: INDEL ###
###################
closeAllConnections()
taxon_split2fasta(vcf_tsv_gain_seq, "INDEL", "Osj(IRGSP)","fasta_annot/indel/gain_indel_Osj_ref.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "INDEL", "Osj","fasta_annot/indel/gain_indel_Osj.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "INDEL","Or","fasta_annot/indel/gain_indel_Or.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "INDEL","Osi","fasta_annot/indel/gain_indel_Osi.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "INDEL","Ob","fasta_annot/indel/gain_indel_Ob.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "INDEL","Og","fasta_annot/indel/gain_indel_Og.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "INDEL","Osj,Osj(IRGSP)","fasta_annot/indel/gain_indel_japonica.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "INDEL","Or,Osj,Osj(IRGSP)","fasta_annot/indel/gain_indel_japonicWild.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "INDEL","Ob,Og","fasta_annot/indel/gain_indel_africa.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "INDEL","Or,Osi,Osj,Osj(IRGSP)","fasta_annot/indel/gain_indel_asia.fasta")

###################
### GAIN: SV ######
###################
closeAllConnections()
taxon_split2fasta(vcf_tsv_gain_seq, "SV", "Osj(IRGSP)","fasta_annot/sv/gain_sv_Osj_ref.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "SV", "Osj","fasta_annot/sv/gain_sv_Osj.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "SV","Or","fasta_annot/sv/gain_sv_Or.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "SV","Osi","fasta_annot/sv/gain_sv_Osi.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "SV","Ob","fasta_annot/sv/gain_sv_Ob.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "SV","Og","fasta_annot/sv/gain_sv_Og.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "SV","Osj,Osj(IRGSP)","fasta_annot/sv/gain_sv_japonica.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "SV","Or,Osj,Osj(IRGSP)","fasta_annot/sv/gain_sv_japonicWild.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "SV","Ob,Og","fasta_annot/sv/gain_sv_africa.fasta")
taxon_split2fasta(vcf_tsv_gain_seq, "SV","Or,Osi,Osj,Osj(IRGSP)","fasta_annot/sv/gain_sv_asia.fasta")

###################
### LOSS: INDEL ###
###################
closeAllConnections()
taxon_split2fasta(vcf_tsv_loss_seq, "INDEL", "Osj(IRGSP)","fasta_annot/indel/loss_indel_Osj_ref.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "INDEL", "Osj","fasta_annot/indel/loss_indel_Osj.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "INDEL","Or","fasta_annot/indel/loss_indel_Or.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "INDEL","Osi","fasta_annot/indel/loss_indel_Osi.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "INDEL","Ob","fasta_annot/indel/loss_indel_Ob.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "INDEL","Og","fasta_annot/indel/loss_indel_Og.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "INDEL","Osj,Osj(IRGSP)","fasta_annot/indel/loss_indel_japonica.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "INDEL","Or,Osj,Osj(IRGSP)","fasta_annot/indel/loss_indel_japonicWild.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "INDEL","Ob,Og","fasta_annot/indel/loss_indel_africa.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "INDEL","Or,Osi,Osj,Osj(IRGSP)","fasta_annot/indel/loss_indel_asia.fasta")

################
### LOSS: SV ###
################
closeAllConnections()
taxon_split2fasta(vcf_tsv_loss_seq, "SV", "Osj(IRGSP)","fasta_annot/sv/loss_sv_Osj_ref.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "SV", "Osj","fasta_annot/sv/loss_sv_Osj.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "SV","Or","fasta_annot/sv/loss_sv_Or.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "SV","Osi","fasta_annot/sv/loss_sv_Osi.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "SV","Ob","fasta_annot/sv/loss_sv_Ob.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "SV","Og","fasta_annot/sv/loss_sv_Og.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "SV","Osj,Osj(IRGSP)","fasta_annot/sv/loss_sv_japonica.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "SV","Or,Osj,Osj(IRGSP)","fasta_annot/sv/loss_sv_japonicWild.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "SV","Ob,Og","fasta_annot/sv/loss_sv_africa.fasta")
taxon_split2fasta(vcf_tsv_loss_seq, "SV","Or,Osi,Osj,Osj(IRGSP)","fasta_annot/sv/loss_sv_asia.fasta")



