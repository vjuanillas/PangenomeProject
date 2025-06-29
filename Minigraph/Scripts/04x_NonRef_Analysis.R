library(dplyr)
library(tidyverse)
# Non_reference analysis
## Quantify the presence of each accession on the non-reference segments: count and lengths
non_reference_counts <- function(pa_matrix) {
  startidx <- which(colnames(pa_matrix) == "SegRank") + 1
  endidx <- ncol(pa_matrix) - 2
  asms <- colnames(pa_matrix)[startidx:endidx]
  
  # initialize dataframe
  non_ref <- pa_matrix[,c("SegmentId", "SegmentLen", "SegRank", asms)]
  
  # filter the ranks
  non_ref <- non_ref %>% filter(SegRank > 0) %>% 
    gather(key = "Accessions", value="Presence", -SegmentId, -SegmentLen, -SegRank) %>%
    filter(Presence > 0)
  
  # to validate, calculate the number of segments where each assembly is present
  non_ref_stats <- non_ref %>% group_by(Accessions)  %>% summarize(count=n(),
                                                                   length=sum(SegmentLen))
return(non_ref_stats)
}

#test1 <- non_reference_counts(AllSegments)
class(AllSegments)

non_reference_sharing <- function(pa_matrix){
  startidx <- which(colnames(pa_matrix) == "SegRank") + 1
  endidx <- ncol(pa_matrix) - 2
  asms <- colnames(pa_matrix)[startidx:endidx]
  
  # initialize data frame
  non_ref <- pa_matrix[,c("SegmentId", "SegmentLen", "SegRank", asms[-1])]
  
  # filter the ranks
  non_ref <- non_ref[non_ref$SegRank>0,]
  
  # update accession names
  # update_asms <-asms[-1] %>% as.character()
  # find combinations per segment
  non_ref$accessions <- apply(non_ref %>% select(asms), 1, 
                              function(x) paste(asms[x > 0], collapse = ","))
  
  # create in long format the non-reference sharing pattern
  non_ref_pattern <- non_ref %>% group_by(accessions) %>% summarize(count = n(), length = sum(SegmentLen))
  
  # ensure that columns are numeric
  non_ref_pattern$count <- as.numeric(non_ref_pattern$count)
  non_ref_pattern$length <- as.numeric(non_ref_pattern$length)
  
  # plot lengths
  plt1 <- ggplot(non_ref_pattern, aes(x = fct_reorder(accessions, length), y = length / 10^6)) +
    geom_col(fill = "blue", color = "black", size = 0.3) + 
    geom_text(aes(label = round(length / 10^6, 2)), hjust = -0.2, size = 2) +  # Adjusted text size for better readability
    scale_y_continuous(limits = c(0, max(non_ref_pattern$length)/10^6 + 5), expand = expansion(mult = c(0, 0.05))) +  # Improved y-axis scaling
    coord_flip() +  # Flip the coordinates for horizontal bar chart
    labs(
      title = "Shared Non-Reference Segment Lengths by Rice accessions",  # Added plot title for context
      x = "Rice Accessions",  # Capitalized for consistency
      y = "Shared Non-Reference Segment Lengths (Mb)"  # Capitalized for consistency
    ) +
    theme_bw(base_size = 9) +  # Use a minimal theme for a cleaner look, adjust base font size for readability
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # Rotate x-axis text for better readability
      axis.text.y = element_text(size = 8),  # Adjust y-axis text size
      plot.title = element_text(hjust = 0.5)  # Center the plot title
    )
  
  plt2 <- ggplot(non_ref_pattern, aes(x = fct_reorder(accessions, count), y = count)) +
    geom_col(fill = "blue", color = "black", size = 0.3) + 
    geom_text(aes(label = count), hjust = -0.2, size = 2) +  # Adjusted text size for better readability
    scale_y_continuous(limits = c(0, max(non_ref_pattern$count) + 5000)) +  # Improved y-axis scaling
    coord_flip() +  # Flip the coordinates for horizontal bar chart
    labs(
      title = "Shared Non-Reference Segment count by rice accessions",  # Added plot title for context
      x = "Rice Accessions",  # Capitalized for consistency
      y = "Shared Non-Reference Segment count"  # Capitalized for consistency
    ) +
    theme_bw(base_size = 9) +  # Use a minimal theme for a cleaner look, adjust base font size for readability
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # Rotate x-axis text for better readability
      axis.text.y = element_text(size = 8),  # Adjust y-axis text size
      plot.title = element_text(hjust = 0.5)  # Center the plot title
    )
  
  print(plt1)
  print(plt2)
}

