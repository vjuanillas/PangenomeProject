#install.packages("gtable", dependencies=TRUE)
#install.packages("ggrepel")
#install.packages("BiocManager")
#install.packages("dendextend")

library(RColorBrewer)
library(ggrepel)
library(tidyverse)
library(dendextend) # for tanglegram
library(ape) # for phylogenetic analysis
library(BiocManager)
library(ggtree) # for phylogenetic tree plotting and annotation
library(gtable)
library(data.table)

# set wd
setwd("/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p90_v2")
# fetch gain and loss files
summed_gain_results <- fread("whole_genome/GAIN_statistics_final.txt")
summed_loss_results <- fread("whole_genome/LOSS_statistics_final.txt")

if (!dir.exists("tree_parsimony")) {
  dir.create("tree_parsimony")
}

##################################
#### Prepare Phylogenetic tree ###
##################################

# data
path_dist_tsv <- "/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p90/all_rice_O_mer.fa.gz.c96f508.community.0.fa.c96f508.eb0f3d3.28057a4.smooth.final.og.dist.tsv"
jaccard_dist_df <- read_tsv(path_dist_tsv) %>%
  arrange(group.a, group.b) %>%
  select(group.a, group.b, jaccard.distance) %>%
  pivot_wider(names_from = group.b, values_from = jaccard.distance) %>%
  column_to_rownames(var = "group.a")

tree <- nj(as.dist(jaccard_dist_df)) # compute the tree using the nj function from the ape package

# Example: Create a dataframe with values mapped to tree nodes
#ggtree(tree) # plot the tree using the ggtree function from the ggtree package
p <- ggtree(
  tree,
  ladderize=T,
  # branch.length
  branch.length="branch.length",
  #layout="daylight"
) + 
  #geom_point(aes(shape=isTip, color=isTip), size=3) + 
  #geom_nodepoint(color="#b5e521", alpha=1/4, size=6) + 
  #geom_tippoint(color="#FDAC4F", shape=8, size=3) +
  geom_tiplab(size=3, color="black") +
  geom_rootedge(rootedge = 0.5) + ggplot2::xlim(-0.1, 0.8) + # to adjust the tree
  geom_treescale(x = 0, y = -1, width = 0.1, fontsize = 3) +  # Add scale bar
  theme_tree2()

print(p)
#################################
### Assign gains and losses #####
#################################
# Merge the dataframes based on CHROM, Taxon, and Type
merged_values <- merge(summed_gain_results, summed_loss_results, by = c("CHROM", "Taxon", "Type"), all = TRUE)
# whole-genome
merged_values_WG <- merged_values[, .(
  Total_Gain_Nbubble = sum(Total_Gain_Nbubble),
  Total_Gain_Length = sum(Total_Gain_Length),
  Total_Loss_Nbubble = sum(Total_Loss_Nbubble),
  Total_Loss_Length = sum(Total_Loss_Length))
  ,by = .(Taxon,Type)] 
merged_values_WG$CHROM <- "WG"

# Add the new column 'coord' based on the analogy
merged_values$node <- with(merged_values, ifelse(Taxon == "Ob", 2, ifelse(Taxon == "Og", 3,ifelse(Taxon == "Or", 4,
                                             ifelse(Taxon == "Osi", 5, ifelse(Taxon == "Osj", 6, ifelse(Taxon == "Osj(IRGSP)", 7, ifelse(Taxon == "Ob,Og", 12,
                                                                         ifelse(Taxon == "Or,Osi,Osj,Osj(IRGSP)", 9, ifelse(Taxon == "Or,Osj,Osj(IRGSP)", 10,
                                                                                       ifelse(Taxon == "Osj,Osj(IRGSP)", 11, NA)))))))))))
merged_values_WG$node <- with(merged_values_WG, ifelse(Taxon == "Ob", 2, ifelse(Taxon == "Og", 3,ifelse(Taxon == "Or", 4,
                                                                                                     ifelse(Taxon == "Osi", 5, ifelse(Taxon == "Osj", 6, ifelse(Taxon == "Osj(IRGSP)", 7, ifelse(Taxon == "Ob,Og", 12,
                                                                                                                                                                                                 ifelse(Taxon == "Or,Osi,Osj,Osj(IRGSP)", 9, ifelse(Taxon == "Or,Osj,Osj(IRGSP)", 10,
                                                                                                                                                                                                                                                    ifelse(Taxon == "Osj,Osj(IRGSP)", 11, NA)))))))))))
fwrite(merged_values_WG, "whole_genome/GAIN_LOSS_statistics_WG_FINAL.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# Plot the tree with node labels
#ggtree(tree, ladderize=TRUE) +
#  geom_text(aes(label=node), hjust=-0.5, size=3) + 
#  theme_tree2()

####################################
### Visualize gains and losses #####
####################################

plot_gain_loss <- function(df, type, chrom, scale, decimal=2, output) {
  node_data <- df %>% filter(Type == type, CHROM == chrom) %>%
    select(node,Taxon,Total_Gain_Nbubble, Total_Gain_Length,Total_Loss_Nbubble,Total_Loss_Length) %>% arrange(node)
  node_data <- node_data %>% mutate(Total_Gain_Length = round(Total_Gain_Length/scale,decimal)) %>% mutate(Total_Loss_Length = round(Total_Loss_Length/scale,decimal)) %>%
    mutate(Gain_Label = paste(Total_Gain_Length, "(",Total_Gain_Nbubble,")"),sep="") %>% mutate(Loss_Label = paste(Total_Loss_Length, "(",Total_Loss_Nbubble,")"))
  # extract node positions in tree
  node_positions <- p$data %>% select(node, x, y) %>% distinct()
  
  # Merge data with node positions
  node_data <- left_join(node_data,node_positions, by="node")
  
  # adjust texts
  external_nodes <- c(1,2,3,4,5,6,7)
  
  viz_p <- p + 
    geom_text(data=node_data, aes(x=x, y=y, label=Gain_Label), 
              nudge_x = ifelse(node_data$node == 7, 0.12, ifelse(node_data$node==11, -0.05,
                                                                 ifelse(node_data$node %in% external_nodes, 0.07, -0.04))), 
              nudge_y = ifelse(node_data$node %in% external_nodes, 0.02, 0.12),
              color="blue", size=3, vjust=-1) +
    geom_text(data=node_data, aes(x=x, y=y, label=Loss_Label), 
              nudge_x = ifelse(node_data$node == 7, 0.12, ifelse(node_data$node==11,-0.05,
                                                                 ifelse(node_data$node %in% external_nodes, 0.07, -0.04))),
              nudge_y = ifelse(node_data$node %in% external_nodes, 0.02, 0.12),
              color="red", size=3, vjust=1)
  
  # Save as PDF
  ggsave(
    filename = paste0(output,"_",type,".png"),
    plot = viz_p,
    width = 10,
    height = 8,
    device = "png",
    dpi = 300
  )
}

for (chr in unique(merged_values$CHROM)){
  plot_gain_loss(merged_values, "SV",chr, 1e6, 2, paste0("tree_parsimony/",chr))
  plot_gain_loss(merged_values, "INDEL", chr, 1e3, 2, paste0("tree_parsimony/",chr))
}

### for whole genome
for (type in unique(merged_values_WG$Type)){
  plot_gain_loss(merged_values_WG, type, "WG", ifelse(type=="SV",1e6,1e3), 2, paste0("tree_parsimony/WG"))
}
