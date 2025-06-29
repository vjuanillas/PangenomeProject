#!/usr/bin/env Rscript
# install.packages("ComplexUpset")
library("ggplot2", warn.conflicts = FALSE)
library("dplyr", warn.conflicts = FALSE)
library("tidyr", warn.conflicts = FALSE)
library("forcats", warn.conflicts = FALSE)
library("magrittr", warn.conflicts = FALSE)
library("UpSetR",warn.conflicts = FALSE)
source("../../Scripts/03_Ref_Analysis.R")

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

nonref_sharing <- function(matrix,decimals,refname) {
  datmat <- read.table(matrix, header = TRUE, stringsAsFactors = FALSE)
  # accessions
  startidx <- which(colnames(datmat) == "SegRank") + 1 # start of name
  endidx <- length(colnames(datmat)) - 2 # end name 
  asms <- colnames(datmat)[startidx:endidx]
  # filter all non-reference without reference column
  datsel <- datmat %>% filter(SegRank > 0) %>% select(-all_of(refname))
  
  datsel$nonrefcol <- apply(datsel %>% select(-StartChrom, -Offset, -total, -accessions,
                                              -SegmentId, -SegmentLen,-SegRank),1, 
                            function(x) paste(asms[x > 0], collapse = ","))
  
  datsum <- datsel %>% group_by(nonrefcol) %>% summarize(shared_count=n(), shared_len=sum(SegmentLen))
  #datsum <- datsel %>% group_by(accessions) %>% summarize(shared_count=n(), shared_len=sum(SegmentLen))
  
  # plot
  theme_set(theme_bw(base_size = 18)+
              theme(panel.grid = element_blank()))
  
  pl1 <- ggplot(datsum,aes(x=fct_reorder(accessions,shared_len),y=shared_len/10^6))+
    geom_col(fill="#56B4E9",col="black",size=0.5)+
    geom_text(aes(label=round(shared_len/10^6,decimals)),hjust=-0.2,size=5)+
    scale_y_continuous(limits=c(0,max(datsum$shared_len)/10^6+5))+
    coord_flip()+
    labs(x="Rice accessions",
         y="Shared non-ref segment lengths (Mb)")
  
  pl2 <- ggplot(datsum,aes(x=fct_reorder(accessions,shared_count),y=shared_count))+
    geom_col(fill="#009E73",col="black",size=0.5)+
    geom_text(aes(label=shared_count),hjust=-0.2,size=5)+
    scale_y_continuous(limits=c(0,max(datsum$shared_count)+5000))+
    coord_flip()+
    labs(x="Rice accessions",
         y="Shared non-ref segment counts")
  
  
  ggsave(pl1, filename = file.path(".", paste0("02_nonref_shared_len.png")),width=12, height=12)
  ggsave(pl2, filename = file.path(".", paste0("02_nonref_shared_count.png")),width=12, height=12)
  
  datsel_inter <- intersection(datsel)
  png("02_nonref_sharing_intersection.png", width = 14, height = 12, units = "in", res = 300)
  theme_set(theme_bw(base_size = 11)+
              theme(panel.grid = element_blank()))
  upset_plot <- upset(fromList(datsel_inter), nsets=6, nintersects = 30, order.by = c("freq"), 
                      set_size.show = FALSE,  mainbar.y.label = "Shared nodes", sets.x.label = "Total nonreference nodes")
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
  nonref_sharing(matrix_file, decimals,refname)
  
  cat("Analysis complete. Plots saved in the current directory.\n")
}

# Run the main function
main()

refname <- "IRGSP"

startidx <- which(colnames(datmat) == "SegRank") + 1 # start of name
endidx <- length(colnames(datmat)) - 2 # end name 
asms <- colnames(datmat)[startidx:endidx][-1]
  # filter all non-reference without reference column
datsel2 <- datmat %>% filter(SegRank > 0) %>% select(-all_of(refname))
  
datsel2$nonrefcol <- apply(datsel2 %>% select(-StartChrom, -Offset, -total, -accessions,
                                              -SegmentId, -SegmentLen,-SegRank),1, 
                            function(x) paste(asms[x > 0], collapse = ","))
  
  datsum <- datsel %>% group_by(nonrefcol) %>% summarize(shared_count=n(), shared_len=sum(SegmentLen))
  #datsum <- datsel %>% group_by(accessions) %>% summarize(shared_count=n(), shared_len=sum(SegmentLen))
  
  # plot
  theme_set(theme_bw(base_size = 18)+
              theme(panel.grid = element_blank()))
  
  pl1 <- ggplot(datsum,aes(x=fct_reorder(accessions,shared_len),y=shared_len/10^6))+
    geom_col(fill="#56B4E9",col="black",size=0.5)+
    geom_text(aes(label=round(shared_len/10^6,decimals)),hjust=-0.2,size=5)+
    scale_y_continuous(limits=c(0,max(datsum$shared_len)/10^6+5))+
    coord_flip()+
    labs(x="Rice accessions",
         y="Shared non-ref segment lengths (Mb)")
  print(pl1)