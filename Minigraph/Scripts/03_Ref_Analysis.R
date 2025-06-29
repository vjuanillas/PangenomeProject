# install.packages("UpSetR")
library(UpSetR)
library(tidyverse)
library(data.table)
library(ggplot2)
library(dplyr)
library(forcats)
# input is still expected to be already read as a dataframe
# All on reference analysis
## determine core, dispensable, and private plus their lengths
reference_categorization <- function(pa_matrix) {
  AllSegments <- pa_matrix
  n_assemblies <- max(AllSegments$total)
  
  # core (present in all assemblies including reference)
  core <- AllSegments %>% filter(total == n_assemblies)
  dispensable <- AllSegments %>% filter(total < n_assemblies & total > 1)
  private <- AllSegments %>% filter(total == 1)
  
  ## counting and adding lengths 
  category_counts <- c(nrow(core), nrow(dispensable), nrow(private), nrow(AllSegments)) 
  category_lengths <- c(sum(core$SegmentLen), sum(dispensable$SegmentLen), sum(private$SegmentLen)
                        , sum(AllSegments$SegmentLen))
  ## Summarize to dataframe and export
  category_stats <- data.frame(category=c("core","dispensable","private","total"), count = category_counts,
                               length = category_lengths)
  write.table(category_stats,file="03_Ref_Category_CountLen.tsv", quote = FALSE,row.names = FALSE)
  write.table(core,file="03_Ref_Category_Core.tsv", quote = FALSE,row.names = FALSE)
  write.table(dispensable,file="03_Ref_Category_Dispensable.tsv", quote = FALSE,row.names = FALSE)
  write.table(private,file="03_Ref_Category_Private.tsv", quote = FALSE,row.names = FALSE)
}

# create an intersection plot
## input: pa_matrix = AllSegments from 02_Segment_matrices.R
## outputs segment_list which can be used to create the plot
intersection <- function(pa_matrix){
  # create a vector of assembly names
  startidx <- which(colnames(pa_matrix) == "SegRank") + 1 # start of name
  endidx <- length(colnames(pa_matrix)) - 2 # end name 
  asms <- colnames(pa_matrix)[startidx:endidx] # slice the dataframe column names
  
  segment_list <- list() # initiate list for intersection
  for (accession in asms) {
    filter_by_segments <- pa_matrix %>% filter(get(accession) != 0 ) %>% pull(SegmentId)
    segment_list[[accession]] <- filter_by_segments 
  }
  
  return(segment_list)
}

# segment quantification of shared nodes
## input: All Segments from 02 
## output: a plot to quantify the segment sharing lengths including the reference 
reference_sharing <- function(pa_matrix) {
  # create a vector of assembly names
  startidx <- which(colnames(pa_matrix) == "SegRank") + 1 
  endidx <- length(colnames(pa_matrix)) -2 
  asms <- colnames(pa_matrix)[startidx:endidx] # slice the data frame column names
  
  sharing_sum <- pa_matrix %>% group_by(accessions) %>% summarize(count = n(),
                                                                  length = sum(SegmentLen))
  sharing_sum$count <- as.numeric(sharing_sum$count)
  sharing_sum$length <- as.numeric(sharing_sum$length)
  
  # Improved plot
  #theme_set(theme_bw(base_size = 9)+
  #            theme(panel.grid = element_blank()))
  
  plt <- ggplot(sharing_sum, aes(x = fct_reorder(accessions, length), y = length / 10^6)) +
    geom_col(fill = "blue", color = "black", size = 0.3) + 
    geom_text(aes(label = round(length / 10^6, 2)), hjust = -0.2, size = 2) +  # Adjusted text size for better readability
    scale_y_continuous(limits = c(0, max(sharing_sum$length)/10^6 + 5), expand = expansion(mult = c(0, 0.05))) +  # Improved y-axis scaling
    coord_flip() +  # Flip the coordinates for horizontal bar chart
    labs(
      title = "Shared Segment Length by Rice accessions including reference",  # Added plot title for context
      x = "Rice Accessions",  # Capitalized for consistency
      y = "Cumulative Segment Length (Mb)"  # Capitalized for consistency
    ) +
    theme_bw(base_size = 9) +  # Use a minimal theme for a cleaner look, adjust base font size for readability
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # Rotate x-axis text for better readability
      axis.text.y = element_text(size = 8),  # Adjust y-axis text size
      plot.title = element_text(hjust = 0.5)  # Center the plot title
    )
  
  return(plt)
}