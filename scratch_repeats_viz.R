library(data.table)
library(tidyverse)
library(dplyr)

setwd("/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p90_v2")
repeats <- fread("whole_genome/repeat_statistics.txt")

repeats$Unmasked_gain <- repeats$Total_Gain_Length - repeats$Masked_Gain_Length 
repeats$Unmasked_loss <- repeats$Total_Loss_Length - repeats$Masked_Loss_Length

repeats_INDEL <- repeats %>% filter(Type == "INDEL")
repeats_SV <- repeats %>% filter(Type == "SV")

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

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
