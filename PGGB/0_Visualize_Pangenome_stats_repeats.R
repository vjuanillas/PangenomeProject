library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

# ===================== #
# REPEATS VISUALIZATION #
# ===================== #
# This script is used to visualize the repeat statistics from the pangenome analysis.
# The script will create a stacked bar plot for the repeat types and their lengths
# across different pangenome categories (Core, Dispensable, Private).
# The script will also create a combined stacked bar plot for the number and length
# of nodes in the pangenome categories.
# The script will save the plots as JPEG files.
# The script assumes that the input data is in the form of a data frame with the following columns:
# - Category: The pangenome category (Core, Dispensable, Private)
# - TE_Type: The type of transposable element (SINEs, LINEs, etc.)
# process the repeats from repeat masker and save as shown in pangenome_statistics.txt

# set working directory
setwd("/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p90")
# load data
repeats <- fread("pangenome_repeats.tsv")
stats <- fread("pangenome_statistics.tsv")

# ============================= #
# Create Visualization function #
# ============================= #

plot_repeats_stacked_bar <- function(data, scale_factor = 1e3, scale_label, output) {
  # Reshape data to long format
  data_long <- data %>%
    pivot_longer(cols = -Category, names_to = "TE_Type", values_to = "Length") %>%
    mutate(Length = Length / scale_factor)  # Scale values
  
  # Define colors, ensuring "Unmasked" is gray
  color_palette <- c(
    "SINEs" = "#E41A1C", "LINEs" = "#377EB8", "Ty1/Copia(LTR)" = "#4DAF4A", 
    "Gypsy/DIRS1(LTR)" = "#984EA3", "DNA Transposons" = "#FF7F00", 
    "Unclassified" = "#FFFF33", "Small RNA" = "#A65628", "Satellites" = "#F781BF", 
    "Simple Repeats" = "magenta", "Low Complexity" = "#66C2A5", 
    "Unmasked" = "gray"
  )
  
  # Create vertical stacked bar plot
  plot <- ggplot(data_long, aes(x = Category, y = Length, fill = TE_Type)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = color_palette) +
    labs(x = "Pangenome Category", y = paste("Length (", scale_label, ")"), fill = "TE Type") +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip() +  # Flip coordinates to make it horizontal
    theme_classic2() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  # Ensure readable labels
  
  # Save as PNG
  ggsave(
    filename = paste0(output, ".jpeg"),
    plot = plot,
    width = 8,
    height = 6 ,
    device = "jpeg",
    dpi = 300
  )
}

plot_stacked_bar_combined <- function(data, output) {
  # Reshape data for stacking (both Number and Length)
  data_long_number <- data %>%
    pivot_longer(cols = "Number", names_to = "Metric", values_to = "Value")
  
  data_long_length <- data %>%
    pivot_longer(cols = "Length (bp)", names_to = "Metric", values_to = "Value")
  
  # Stacked Bar Plot: Based on Number
  p1 <- ggplot(data_long_number, aes(y = Categories, x = Value, fill = Categories)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = c("Core" = "#4DAF4A", "Dispensable" = "#377EB8", "Private" = "#E41A1C")) +
    labs(x = "Count", y = NULL, title = "Number of Nodes") +
    theme_classic()
  
  # Stacked Bar Plot: Based on Length
  p2 <- ggplot(data_long_length, aes(y = Categories, x = Value/1e6, fill = Categories)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = c("Core" = "#4DAF4A", "Dispensable" = "#377EB8", "Private" = "#E41A1C")) +
    labs(x = "Length (Mb)", y = NULL, title = "Length") +
    theme_classic()
  
  # Arrange both plots together
  combined_plot <- ggarrange(p1, p2, ncol = 1, common.legend = TRUE, legend = "right")
  
  # Save as JPEG
  ggsave(filename = paste0(output, ".jpeg"), 
         plot = combined_plot, 
         width = 8, 
         height = 10, 
         dpi = 300)
}

# plot the stacked bar plot for repeats
plot_repeats_stacked_bar(repeats, 1e6, "Mb", "pangenome_repeats")
plot_stacked_bar_combined(stats,"pangenome_stats")
