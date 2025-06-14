!/usr/bin/env Rscript

library(tidyverse)
library(dendextend) # for tanglegram
library(ape) # for phylogenetic analysis
library(ggtree) # for phylogenetic tree plotting and annotation

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if we have exactly 3 arguments
if(length(args) != 3) {
  stop("Usage: Rscript make_tree.R path_dist_tsv plot_title output_file_path", call. = FALSE)
}

# Assign arguments
path_dist_tsv <- args[1]
plot_title <- args[2]
output_file_path <- args[3]

# Read matrices
jaccard_dist_df <- read_tsv(path_dist_tsv) %>%
  arrange(group.a, group.b) %>%
  select(group.a, group.b, jaccard.distance) %>%
  pivot_wider(names_from = group.b, values_from = jaccard.distance) %>%
  column_to_rownames(var = "group.a")

#euclidean_dist_df <- read_tsv(path_dist_tsv) %>%
#  arrange(group.a, group.b) %>%
#  select(group.a, group.b, euclidean.distance) %>%
#  pivot_wider(names_from = group.b, values_from = euclidean.distance) %>%
#  column_to_rownames(var = "group.a")

jaccard_hc <- as.dist(jaccard_dist_df) %>% hclust()
#euclidean_hc <- as.dist(euclidean_dist_df) %>% hclust()

#euclidean_hc <- as.dist(euclidean_dist_df) %>% hclust()

# Plot
pdf(output_file_path) # Start PDF device, replace with svg(output_file_path) for SVG format
plot(
  jaccard_hc,
  # Label at same height
  hang = -1,
  main = plot_title,
  xlab = 'Haplotype',
  ylab = 'Jaccard distance',
  sub = '',
  cex = 1.8,       # Adjusts the size of points/text in the plot
  cex.lab = 1.6,   # Adjusts the size of x and y labels
  cex.axis = 1.6,  # Adjusts the size of axis text
  cex.main = 1.6,  # Adjusts the size of the main title
  cex.sub = 1.6,    # Adjusts the size of the subtitle
  lwd = 2  # Increases the width of the branch lines
)
dev.off() # Close the device



library(RColorBrewer)
library(ggrepel)

tree <- nj(as.dist(jaccard_dist_df)) # compute the tree using the nj function from the ape package

# Plot and save the NJ tree using ggtree
output_file_nj <- gsub(".pdf$", "_nj_tree.pdf", output_file_path)  # Modify output file name for NJ tree

pdf(output_file_nj)  # Save as PDF
#ggtree(tree) # plot the tree using the ggtree function from the ggtree package
ggtree(
  tree,
  ladderize=T,

  # branch.length
  branch.length="branch.length",

  #layout="daylight"
) +
  #geom_point(aes(shape=isTip, color=isTip), size=3) + 
  geom_nodepoint(color="#b5e521", alpha=1/4, size=6) +
  #geom_tippoint(color="#FDAC4F", shape=8, size=3) +
  geom_tiplab(size=3, color="black") +
  geom_rootedge(rootedge = 0.5) + ggplot2::xlim(0, 0.8) +
  geom_treescale(x = 0, y = -1, width = 0.1, fontsize = 3) +  # Add scale bar
  theme_tree2()

dev.off()