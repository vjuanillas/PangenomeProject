#install.packages("ggpattern")
# install.packages("ggpubr")
library(ggplot2)
library(ggpattern)
library(dplyr)
library(readr)
library(stringr)  # For regex-based text processing
library(tidyr)
library(ggpubr)

setwd("/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p90_v3/whole_genome")

###
## # add data value in the middle of the donuts
# operate from here
# Read gains and loss data with SV and INDELs
gain_raw <- fread("Gain_parsimony_all_summary.txt")
loss_raw <- fread("Loss_parsimony_all_summary.txt")


#######################################################################
## visualize distribution of parsimonious and non-parsimonious cases ##
#######################################################################

# Summarize data by subtype, par_not status
gain_summary <- gain_raw[, .(Count = sum(N_bubble), Length = sum(GAIN_length), Ave_length = sum(GAIN_length)/sum(N_bubble)), 
                         by = .(CHROM, subtype, par_not)]
loss_summary <- loss_raw[, .(Count = sum(N_bubble), Length = sum(LOSS_length), Ave_length = sum(LOSS_length)/sum(N_bubble)),
                         by = .(CHROM, subtype, par_not)]

# Assign mutation type and category
ins_indel <- gain_summary %>% filter(subtype == "INS") %>% mutate(type = "INS", Category = "INDEL")
ins_sv <- gain_summary %>% filter(subtype == "SV_INS") %>% mutate(type = "INS", Category = "SV")
del_indel <- loss_summary %>% filter(subtype == "DEL") %>% mutate(type = "DEL", Category = "INDEL")
del_sv <- loss_summary %>% filter(subtype == "SV_DEL") %>% mutate(type = "DEL", Category = "SV")

# Combine datasets
all_data <- bind_rows(del_indel, ins_indel, del_sv, ins_sv)

# Create a combined factor for type and parsimony
all_data$type_parsimony <- paste(all_data$type, all_data$par_not, sep = "_")

# Remove 'chr' prefix and convert to factor
all_data$CHROM <- str_remove(all_data$CHROM, "^chr")
all_data$CHROM <- factor(all_data$CHROM, levels = 1:12, ordered = TRUE)

# Create a combined factor for type and parsimony
WG_var_stats <- all_data[,.(WG_count = sum(Count), WG_length = sum(Length)),
                 by=.(subtype, type_parsimony)]
# Compute total WG_count per subtype
WG_var_stats[, Total_count := sum(WG_count), by = subtype]

# Compute percentage for each type_parsimony within its subtype
WG_var_stats[, Percent := (WG_count / Total_count) * 100]
WG_var_stats$Percent <- round(WG_var_stats$Percent, 2)

# Create a custom color palette that varies by both type and parsimony
fill_colors <- c(
  "DEL_par" = "#B22222",    # Dark red for parsimonious DEL
  "DEL_not" = "#f7b5e9",    # Light pink for non-parsimonious DEL
  "INS_par" = "#1E3A8A",    # Dark blue for parsimonious INS
  "INS_not" = "#ADD8E6"     # Light blue for non-parsimonious INS
)

# add data value in the middle
donut_par_not <- ggdonutchart(WG_var_stats, "Percent", fill = "type_parsimony", label = "type_parsimony") +
  facet_wrap(~ subtype) +  # Separate plots for GAINS and LOSS
  scale_fill_manual(values = fill_colors,
                    labels = c("DEL (Non-parsimonious)", "DEL (Parsimonious)", 
                               "INS (Non-parsimonious)", "INS (Parsimonious)")) + 
  labs(fill = "Type") +
  theme_void() +
  theme(legend.position = "right")

# Save as PDF
ggsave(
  filename = paste0("parsi_combination_distribution.jpeg"),
  plot = donut_par_not,
  width = 10,
  height = 6,
  device = "jpeg",
  dpi = 300
)
################################
### Bar plots per chromosome ###
################################

# Plot stacked bar chart for number
plot_count <- ggplot(all_data, aes(x = CHROM, y = Count, fill = type_parsimony)) +
  geom_bar(stat = "identity", position = "stack", colour = "black", linewidth = 0.3) +
  facet_wrap(~ Category, scales = "free") +
  scale_fill_manual(values = fill_colors,
                    labels = c("DEL (Non-parsimonious)", "DEL (Parsimonious)", 
                               "INS (Non-parsimonious)", "INS (Parsimonious)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Chromosome", 
       y = "Number of Variants", 
       fill = "Type") +
  theme_classic() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 12, face = "bold"))

# Plot stacked bar chart for total length
plot_totalLength <- ggplot(all_data, aes(x = CHROM, y = Length/1e6, fill = type_parsimony)) +
  geom_bar(stat = "identity", position = "stack", colour = "black", linewidth = 0.3) +
  facet_wrap(~ Category, scales = "free") +
  scale_fill_manual(values = fill_colors,
                    labels = c("DEL (Non-parsimonious)", "DEL (Parsimonious)", 
                               "INS (Non-parsimonious)", "INS (Parsimonious)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Chromosome", 
       y = "Cumulative Length (Mb)", 
       fill = "Type") +
  theme_classic() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 12, face = "bold"))

# stacked bar chart for ave length
plot_avgLength <- ggplot(all_data, aes(x = CHROM, y = Ave_length, fill = type_parsimony)) +
  geom_bar(stat = "identity", position = "stack", colour = "black", linewidth = 0.3) +
  facet_wrap(~ Category, scales = "free") +
  scale_fill_manual(values = fill_colors,
                    labels = c("DEL (Non-parsimonious)", "DEL (Parsimonious)", 
                               "INS (Non-parsimonious)", "INS (Parsimonious)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Chromosome", 
       y = "Average Length (bp)", 
       fill = "Type") +
  theme_classic() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 12, face = "bold"))

combined_plot <- ggarrange(plot_count, plot_totalLength, plot_avgLength)
# Save as PDF
ggsave(
  filename = paste0("variant_distribution_domRice.pdf"),
  plot = combined_plot,
  width = 36,
  height = 24,
  device = "pdf",
  dpi = 300
)
######################################################
### Split plot of non-parsimonious patterns ###
######################################################

gain_raw_nonpar <- gain_raw %>% filter(par_not == "not")
gain_raw_nonpar <- gain_raw_nonpar[,.(Total_count = sum(N_bubble), Total_length = sum(GAIN_length)),
                by=.(ALT_accn, subtype)]

loss_raw_nonpar <- loss_raw %>% filter(par_not == "not")
loss_raw_nonpar <- loss_raw_nonpar[,.(Total_count = sum(N_bubble), Total_length = sum(LOSS_length)),
                                   by=.(ALT_accn, subtype)]

gain_raw_nonpar_INS <- gain_raw_nonpar %>% filter(subtype == "INS")
gain_raw_nonpar_SVINS <- gain_raw_nonpar %>% filter(subtype == "SV_INS")
loss_raw_nonpar_DEL <- loss_raw_nonpar %>% filter(subtype == "DEL")
loss_raw_nonpar_SVDEL <- loss_raw_nonpar %>% filter(subtype == "SV_DEL")


create_two_sided_barplot <- function(data, plot_name, unit = "Kbp", scale_factor = 1e3, color="#ADD8E6") {
  # Input validation
  required_cols <- c("ALT_accn", "Total_count", "Total_length")
  if (!all(required_cols %in% colnames(data))) {
    stop("Data must contain columns: ALT_accn, Total_count, Total_length")
  }
  
  # Filter for top 10 based on total length
  data_top10 <- data %>%
    arrange(desc(Total_length)) %>%
    slice_head(n = 10)
  
  axis_margin <- 5.5
  
  # Create left plot
  p1 <- ggplot(data_top10, aes(Total_count, ALT_accn)) +
    geom_col(fill = color, colour = "black", linewidth = 0.3) +
    scale_x_reverse() +
    scale_y_discrete(position = "right") +
    labs(x = "Total Count") +
    ggtitle(plot_name) + 
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(axis_margin, 0, axis_margin, axis_margin),
      axis.title.x = element_text(hjust = 0.5),
      plot.title = element_text(vjust = 0.5, face = "bold")
    )
  
  # Create right plot
  p2 <- ggplot(data_top10, aes(Total_length/scale_factor, ALT_accn)) +
    geom_col(fill = color, colour = "black", linewidth = 0.3) +
    labs(x = paste0("Total Length (", unit, ")")) +
    ggtitle(" ") +
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      plot.margin = margin(axis_margin, axis_margin, axis_margin, 0),
      axis.text.y.left = element_text(margin = margin(0, axis_margin, 0, axis_margin),
                                      hjust = 0.5),
      axis.title.x = element_text(hjust = 0.5)
    )
  
  # Combine plots and save as PDF
  combined_plot <- ggarrange(p1, p2)
  
  # Save as PDF
  ggsave(
    filename = paste0(plot_name, ".jpeg"),
    plot = combined_plot,
    width = 12,
    height = 6,
    device = "jpeg",
    dpi = 300
  )
  
  # Return the plot object in case it's needed
  return(combined_plot)
}

create_two_sided_barplot(gain_raw_nonpar_INS, "Non-parsimonious Insertions")
create_two_sided_barplot(gain_raw_nonpar_SVINS, "Non-parsimonious SV Insertions", scale_factor = 1e6, unit = "Mbp")
create_two_sided_barplot(loss_raw_nonpar_DEL, "Non-parsimonious Deletions", color="#f7b5e9")
create_two_sided_barplot(loss_raw_nonpar_SVDEL, "Non-parsimonious SV Deletions", color="#f7b5e9", scale_factor = 1e6, unit = "Mbp")
