#!/usr/bin/env Rscript
# install.packages("ComplexUpset")
library("ggplot2", warn.conflicts = FALSE)
library("dplyr", warn.conflicts = FALSE)
library("tidyr", warn.conflicts = FALSE)
library("forcats", warn.conflicts = FALSE)
library("magrittr", warn.conflicts = FALSE)
library("UpSetR",warn.conflicts = FALSE)
source("/Users/jongpaduhilao/Desktop/LAB_Files/Minigraph/Trial_4/Scripts/03_Ref_Analysis.R")

# takes 02_Presence_Absence_ALL.tsv 
# nonreference sequences are extracted
# returns a dataframe of assembly and non-ref seq
calculate_nonref <- function(matrix,refname) {
  datmat <- read.table(matrix, header = TRUE, stringsAsFactors = FALSE)
  
  # select non-ref nodes
  # convert to long format
  
  datsel <- datmat %>% filter(SegRank >0) %>%
    select(-StartChrom, -Offset, -total, -accessions) %>%  # Remove columns not needed in long format
    gather(key = "accession",
           value = "presence",
           -SegRank, -SegmentLen, -SegmentId) %>% 
    filter(presence > 0)
  
# obtain non-ref count and length
  datout <- datsel %>% group_by(accession) %>%
    summarize(no_segment=n(), 
              length_segment = sum(SegmentLen))
  refassemb  <- refname
  datout   <- datout  %>% filter(accession != refassemb)   
  
  
  #outfile  <- file.path(outbase, paste0(graphtype, "_nonref_analysis.tsv"))
  
  #write.table(datout, file=outfile, quote = FALSE, row.names=FALSE)
  return(datout)
}

nonref_sharing <- function(matrix, decimals) {
  datmat <- read.table(matrix, header = TRUE, stringsAsFactors = FALSE)
  
  # filter all non-reference without reference column
  datsel <- datmat %>% filter(SegRank > 0)
  datsum <- datsel %>%
    group_by(accessions) %>%
    summarize(shared_count = n(), shared_len = sum(SegmentLen)) %>%
    ungroup()
  
  # Keep only top 20 by shared_len
  top20 <- datsum %>% slice_max(shared_len, n = 20)
  
  # plot
  theme_set(theme_bw(base_size = 18) +
              theme(panel.grid = element_blank()))
  
  pl1 <- ggplot(top20, aes(x = fct_reorder(accessions, shared_len), y = shared_len / 10^6)) +
    geom_col(fill = "#56B4E9", col = "black", size = 0.5) +
    geom_text(aes(label = round(shared_len / 10^6, decimals)), hjust = -0.2, size = 5) +
    scale_y_continuous(limits = c(0, max(top20$shared_len) / 10^6 + 5)) +
    coord_flip() +
    labs(x = "Rice accessions",
         y = "Shared non-ref segment lengths (Mb)")
  
  # Keep only top 20 by shared_count
  top20_count <- datsum %>% slice_max(shared_count, n = 20)
  
  pl2 <- ggplot(top20_count, aes(x = fct_reorder(accessions, shared_count), y = shared_count)) +
    geom_col(fill = "#009E73", col = "black", size = 0.5) +
    geom_text(aes(label = shared_count), hjust = -0.2, size = 5) +
    scale_y_continuous(limits = c(0, max(top20_count$shared_count) + 5000)) +
    coord_flip() +
    labs(x = "Rice accessions",
         y = "Shared non-ref segment counts")
  
  ggsave(pl1, filename = file.path(".", paste0("02_nonref_shared_len_top20.png")), width = 12, height = 12)
  ggsave(pl2, filename = file.path(".", paste0("02_nonref_shared_count_top20.png")), width = 12, height = 12)
  
  datsel_inter <- intersection(datsel)
  png("02_nonref_sharing_intersection.png", width = 14, height = 12, units = "in", res = 300)
  theme_set(theme_bw(base_size = 9) +
              theme(panel.grid = element_blank()))
  upset_plot <- upset(fromList(datsel_inter), nsets = 6, nintersects = 35, order.by = c("freq"),
                      set_size.show = FALSE, text.scale = 1.5, mainbar.y.label = "Shared nodes", sets.x.label = "Total nonreference nodes")
  print(upset_plot)
  dev.off()
}

# Main execution
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    stop("Usage: Rscript script.R <matrix_file> [refname] [decimal_places]")
  }
  
  matrix_file <- args[1]
  refname <- if (length(args) > 1) args[2] else "IRGSP"
  decimals <- if (length(args) > 2) as.numeric(args[3]) else 2
  
  cat("Calculating non-reference sequences...\n")
  nonref_data <- calculate_nonref(matrix_file, refname)
  print(head(nonref_data))
  
  cat("Generating plots...\n")
  nonref_sharing(matrix_file, decimals)
  
  cat("Analysis complete. Plots saved in the current directory.\n")
}

# Run the main function
main()