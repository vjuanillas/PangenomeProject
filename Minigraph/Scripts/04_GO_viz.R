library(ggplot2)
library(reshape2)
library(dplyr)

# input: a long dataframe with 6 columns: number, gene ontology terms, 
# number of genes, percent gene hits, percent function hits, and category
# output: a plot for visualization
go_visualize_1 <- function(df) {
  df_mutated <- df %>% mutate(category = factor(category, levels = c("core", "dispensable", "private")))
  plt <- ggplot(df_mutated, aes(x = category, y = `Gene Ontology`)) +
    geom_point(aes(size = `Number of genes`, color = category)) +
    scale_size_continuous(range = c(1, 10)) +
    scale_color_manual(values = c("core" = "darkgreen", "dispensable" = "darkblue", "private" = "darkred")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    labs(x = "Category", y = "Gene Ontology Terms", size = "Number of genes", ) +
    guides(color = FALSE)
return(plt)
}

go_visualize_2 <- function(df) {
  df_mutated <- df %>% mutate(category = factor(category, levels = unique(df$category)))
  plt <- ggplot(df_mutated, aes(x = category, y = `Gene Ontology`)) +
    geom_point(aes(size = `Number of genes`, color = category)) +
    scale_size_continuous(range = c(1, 10)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    labs(x = "Category", y = "Gene Ontology Terms", size = "Number of genes", ) +
    guides(color = FALSE)
return(plt)
}