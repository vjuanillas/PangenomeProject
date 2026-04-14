library(jsonlite)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(gridExtra)

# convert JSON to table
parse_panther_json <- function(file, group_name, ontology) {
  
  json_data <- fromJSON(file, flatten = TRUE)
  
  df <- json_data$results$result %>%
    as.data.frame()
  
  df_clean <- df %>%
    transmute(
      term_id = term.id,
      term_label = term.label,
      fold_enrichment = suppressWarnings(as.numeric(fold_enrichment)),
      p_value = suppressWarnings(as.numeric(pValue)),
      fdr = suppressWarnings(as.numeric(fdr)),
      gene_count = suppressWarnings(as.numeric(number_in_list)),
      group = group_name,
      ontology = ontology
    ) %>%
    
  # Filter out the Unclassified
  filter(
    !is.na(term_label),
    !str_detect(tolower(term_label), "unclassified")
  ) %>%
    
    mutate(
      fdr = ifelse(is.na(fdr) | fdr <= 0, 1e-300, fdr),
      gene_count = ifelse(is.na(gene_count), 0, gene_count)
    )
  
  return(df_clean)
}

# Function to plot pan and variable
plot_bubble <- function(df, groups_to_keep, title, top_n = 10) {
  
  plot_df <- df %>%
    filter(group %in% groups_to_keep) %>%
    group_by(group, ontology) %>%
    arrange(fdr, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(
      term_label = str_wrap(term_label, width = 45),
      log_fdr = -log10(fdr),
      term_label = factor(term_label, levels = rev(unique(term_label)))
    )
  
  ggplot(plot_df,
         aes(x = group, y = term_label)) +
    
    geom_point(aes(size = gene_count,
                   color = log_fdr),
               alpha = 0.85) +
    
    facet_wrap(~ontology, scales = "free_y") +
    
    scale_size_continuous(name = "Gene Count") +
    scale_color_gradient(name = "-log10(FDR)",
                         low = "skyblue",
                         high = "red") +
    
    labs(title = title, x = "", y = "GO Terms") +
    
    theme_minimal(base_size = 16) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 16),
      panel.grid.major.y = element_blank()
    )
}

# Function to plot core, disp, and private
plot_single_group <- function(df, group_name, title, top_n = 10) {
  
  plot_df <- df %>%
    filter(group == group_name) %>%
    group_by(ontology) %>%
    arrange(fdr, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(
      term_label = str_wrap(term_label, width = 45),
      log_fdr = -log10(fdr),
      term_label = factor(term_label, levels = rev(unique(term_label)))
    )
  
  ggplot(plot_df,
         aes(x = group, y = term_label)) +
    
    geom_point(aes(size = gene_count,
                   color = log_fdr),
               alpha = 0.85) +
    
    facet_wrap(~ontology, scales = "free_y") +
    
    scale_size_continuous(name = "Gene Count") +
    scale_color_gradient(name = "-log10(FDR)",
                         low = "skyblue",
                         high = "red") +
    
    labs(title = title, x = "", y = "GO Terms") +
    
    theme_minimal(base_size = 16) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 16),
      panel.grid.major.y = element_blank()
    )
}

# what's the top 10 gene terms?
run_panther_pipeline <- function(file_list,
                                 output_prefix = "panther",
                                 top_n = 10,
                                 keep_ontologies = c("BP")) {
  
  file_list <- file_list[names(file_list) %in% keep_ontologies]
  
  all_results <- list()
  
  for (ont in names(file_list)) {
    
    files <- file_list[[ont]]
    
    core <- parse_panther_json(files$core, "core", ont)
    disp <- parse_panther_json(files$disp, "dispensable", ont)
    priv <- parse_panther_json(files$priv, "private", ont)
    
    combined <- bind_rows(core, disp, priv)
    
    all_results[[ont]] <- combined
    
    write_tsv(combined,
              paste0(output_prefix, "_", ont, ".tsv"))
  }
  
  final_df <- bind_rows(all_results)
  
  write_tsv(final_df,
            paste0(output_prefix, "_ALL.tsv"))
  
  #display individual bubble plots
  plot_core <- plot_single_group(final_df, "core",
                                 "Core Genome GO Enrichment",
                                 top_n)
  
  plot_disp <- plot_single_group(final_df, "dispensable",
                                 "Dispensable Genome GO Enrichment",
                                 top_n)
  
  plot_priv <- plot_single_group(final_df, "private",
                                 "Private Genome GO Enrichment",
                                 top_n)
  
  # plot the pan and variable
  plot_pan <- plot_bubble(final_df,
                          c("core", "dispensable", "private"),
                          "Pan Genome: Core vs Dispensable vs Private",
                          top_n)
  
  plot_variable <- plot_bubble(final_df,
                               c("dispensable", "private"),
                               "Variable Genome: Dispensable vs Private",
                               top_n)
  
  # Save to PDF all of the bubble plots
  pdf(paste0(output_prefix, "_bubble_plots.pdf"),
      width = 14, height = 12)
  
  print(plot_core)
  print(plot_disp)
  print(plot_priv)
  print(plot_pan)
  print(plot_variable)
  
  dev.off()
  
  return(list(
    full_data = final_df,
    core_plot = plot_core,
    dispensable_plot = plot_disp,
    private_plot = plot_priv,
    pan_plot = plot_pan,
    variable_plot = plot_variable
  ))
}


files <- list(
  BP = list(
    core = "core_panther_enrichment_bp.json",
    disp = "dispensable_panther_enrichment_bp.json",
    priv = "private_panther_enrichment_bp.json"
  )
)

result <- run_panther_pipeline(
  files,
  output_prefix = "rice_panther_BP",
  top_n = 10,
  keep_ontologies = c("BP")
)

print(result$core_plot)
print(result$dispensable_plot)
print(result$private_plot)