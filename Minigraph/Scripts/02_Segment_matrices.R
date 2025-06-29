#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(dplyr)

# this script should be able to create a presence-absence matrix of segment usage
args  <- commandArgs(trailingOnly = TRUE)
segmentAll <- fread(args[1], header = TRUE, stringsAsFactors = FALSE)
segmentAll <- as.data.frame(segmentAll)

# assign asm name based on the first input 
# the reason why I decided to assign the assemblies name like this is because to make debugging faster
# an error in terms of names would always end up in the previous script as no name assignment was done here
start_idx <- which(colnames(segmentAll) == "SegRank") + 1
asms <- colnames(segmentAll)[start_idx:length(colnames(segmentAll))]

# create a dataframe for assemblies plus their ranks
asms_rank <- data.frame(rankid=seq(0,length(asms)-1),
                        accessions=asms,stringsAsFactors = FALSE)

# Based on coverage only
# evaluate whether values under each accession number is 0 or 1 and assign a value
segmentAll[, asms] <- lapply(segmentAll[, asms], function(x) ifelse(x>0, 1, 0)) 

# Based on rank
# update the matrix such that an assembly whose rank is equal to the node's rank is certainly tagged as 1
segmentAll_2 <- segmentAll

for (i in seq(1,nrow(segmentAll))){
  rank <- segmentAll[i,"SegRank"]
  accession_id <- asms_rank[asms_rank$rankid==rank,"accessions"]
  segmentAll_2[i,accession_id] <- 1
}

# Add total column by calculating the number of accessions per segment then label what accessions are present under the accession column
segmentAll_2$total <- rowSums(segmentAll_2[,asms])
matrix_only <- segmentAll_2 # to save the matrix only

# put labels to each of the segments
segmentAll_2 <- segmentAll_2 %>% rowwise() %>% 
  mutate(accessions = paste(asms[which(c_across(all_of(asms))== 1)], collapse = ",")) %>%
  ungroup()

names_only <- data.frame(SegmentId=segmentAll_2$SegmentId, Accessions=segmentAll_2$accessions)


# write to csv
##output files: 
### all tsv
write.table(segmentAll_2,file="02_Presence_Absence_ALL.tsv", quote = FALSE,row.names = FALSE)
### matrix only
write.table(matrix_only,file="02_Presence_Absence_Matrix_only.tsv", quote = FALSE,row.names = FALSE)
### names only
write.table(names_only,file="02_Presence_Absence_Names_only.tsv", quote = FALSE,row.names = FALSE)