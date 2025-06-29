library(data.table) # for fread
library(tidyverse)
library(dplyr)

data <- "/Users/jongpaduhilao/Desktop/LAB_Files/Initial_Pangenome_analysis/Trial_4/Raw_data/data/05x4_biallelic_asm2annot.intersect.tsv"
data <- fread(data,sep = "\t", header = FALSE)

data$V5 <- ifelse(data$V6 == 1, "Insertion", data$V5) 
data$V5 <- ifelse(data$V7 == 1, "Deletion", data$V5)

summarized <- data %>% group_by(V5) %>%
  summarise(
    Count = n(),
    CumLengthShort = sum(V6, na.rm = TRUE),
    CumLengthLong = sum(V7, na.rm = TRUE)
  )

sum(summarized$Count)

data %>% filter(V12 == "CDS") %>% count()

# quantifying lengths of variants over gene structures/features
# retrieving matching ID values
data_gene_structure <- "/Users/jongpaduhilao/Desktop/LAB_Files/Initial_Pangenome_analysis/Trial_4/Raw_data/data/06_biallelic_sv_annot.tsv"
data_gene_structure <- fread(data_gene_structure, sep = "\t", header=FALSE)

matched_gene_structures <- data_gene_structure %>% left_join(data %>% distinct(V4, .keep_all = TRUE) %>%
                                                               select(V4,V6, V7), by = c("V1" = "V4" ))
# write.table(matched_gene_structures, file="06_matched_gene_structures_withLengths.tsv", sep="\t", row.names = FALSE, quote=FALSE)

summarized_gene_struc <- matched_gene_structures %>% group_by(V3) %>%
  summarise(
    n_Insertions = sum(V2 == "Insertion"),
    length_insertions = sum(V7[V2 == "Insertion"], na.rm = TRUE),
    n_Deletions = sum(V2 == "Deletion"),
    length_deletion = sum(V6[V2 == "Deletion"], na.rm = TRUE)
  )

matched_gene_structures %>% group_by(V3) %>%
  summarise(
    n_Insertions = sum(V2 == "Insertion"),
    n_Deletions = sum(V2 == "Deletions"),
    n_AltIns = sum(V2 == "AltIns"),
    n_AltDel = sum(V2 == "AltDel")
  )
