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

#####################
### VISUALIZATION ###
#####################

# process the repeats from repeat masker and save as shown in repeat_statistics.txt

library(data.table)
library(tidyverse)
library(dplyr)
library(forcats)

setwd("/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p90_v2")
repeats <- fread("whole_genome/repeat_statistics.txt")

repeats$Unmasked_gain <- repeats$Total_Gain_Length - repeats$Masked_Gain_Length 
repeats$Unmasked_loss <- repeats$Total_Loss_Length - repeats$Masked_Loss_Length
fwrite(repeats, "whole_genome/Repeats_stats.txt", sep = "\t", quote=F, row.names = F, col.names=T)

repeats_INDEL <- repeats %>% filter(Type == "INDEL")
repeats_SV <- repeats %>% filter(Type == "SV")

plot_stacked_bar <- function(data, scale_factor = 1e3,scale_label,output) {
  # Define categories for gain and loss
  gain_cols <- c("SINEs_gain", "LINEs_gain", "Ty1/Copia(LTR)_gain", "Gypsy/DIRS1(LTR)_gain", 
                 "DNA Transposons_gain", "Unclassified_gain", "Small RNA_gain", "Satellites_gain", 
                 "Simple Repeats_gain", "Low Complexity_gain", "Unmasked_gain")
  
  loss_cols <- c("SINEs_loss", "LINEs_loss", "Ty1/Copia(LTR)_loss", "Gypsy/DIRS1(LTR)_loss", 
                 "DNA Transposons_loss", "Unclassified_loss", "Small RNA_loss", "Satellites_loss", 
                 "Simple Repeats_loss", "Low Complexity_loss", "Unmasked_loss")
  
  # Reshape data to long format
  data_long <- data %>% 
    pivot_longer(cols = c(all_of(gain_cols), all_of(loss_cols)), 
                 names_to = "Category", values_to = "Length") %>% 
    mutate(Length = Length / scale_factor, 
           Type = ifelse(grepl("_gain", Category), "Gain", "Loss"),
           Category = gsub("_gain|_loss", "", Category))
  
  # Define colors, making sure Unmasked_Loss_Length is gray
  color_palette <- c(
    "SINEs" = "#E41A1C", "LINEs" = "#377EB8", "Ty1/Copia(LTR)" = "#4DAF4A", 
    "Gypsy/DIRS1(LTR)" = "#984EA3", "DNA Transposons" = "#FF7F00", 
    "Unclassified" = "#FFFF33", "Small RNA" = "#A65628", "Satellites" = "#F781BF", 
    "Simple Repeats" = "magenta", "Low Complexity" = "#66C2A5", 
    "Unmasked" = "gray")
  
  # Ensure category order is maintained
  data_long$Category <- factor(data_long$Category, levels = gsub("_gain|_loss", "", gain_cols))
  
  # Adjust taxon labels to create two bars per taxon
  data_long <- data_long %>% 
    mutate(Taxon = paste(Taxon, Type, sep = "_"))
  
  # Determine color of x-axis labels
  taxon_colors <- ifelse(grepl("_Loss", levels(factor(data_long$Taxon))), "red", "blue")
  
  # Create plot
  plot <- ggplot(data_long, aes(x = Taxon, y = Length, fill = Category)) +
    geom_bar(stat = "identity", position = position_stack(reverse=TRUE),
             color="black",linewidth=0.2) +
    scale_fill_manual(values = color_palette) +
    coord_flip() +
    labs(x = "Taxon", y = paste("Length (", scale_label, ")"), fill = "Category") +
    scale_y_continuous(expand=c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(color=taxon_colors)) 
  # Save as PDF
  ggsave(
    filename = paste0(output,".png"),
    plot = plot,
    width = 10,
    height = 8,
    device = "png",
    dpi = 300
    
  )
}


plot_stacked_bar(repeats_INDEL, 1e3, "Kb", "whole_genome/Repeats_INDEL")
plot_stacked_bar(repeats_SV, 1e6, "Mb", "whole_genome/Repeats_SV")


