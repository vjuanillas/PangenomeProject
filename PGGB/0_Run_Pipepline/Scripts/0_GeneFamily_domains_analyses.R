library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)
library(tidyr)

# import merged uniprot data
setwd("/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p60/fasta_annot/tsv/")
uniprot_gene <- read.csv("/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p60/fasta_annot/tsv/uniprot/all_categories_uniprot.tsv", 
                         sep ="\t",header=TRUE)
colnames(uniprot_gene) <- c("From", "Entry","Reviewed", "Entry Name","Protein names","Gene Names","Organism","Length","Protein families","Domain [FT]",
                            "Sample")


protein_family_function <- function(df, sample) {
  counts <- df %>%
    #distinct(`Protein families`, .keep_all = TRUE) %>%
    # Split the "Protein families" column by commas
    mutate(family_split = str_split(`Protein families`, "\\s*,\\s*|\\s*;\\s*")) %>%
    unnest(family_split) %>%
    # Clean whitespace
    mutate(family_split = str_trim(family_split)) %>%
    # Filter out anything containing "superfamily" or "subfamily" (case insensitive)
    filter(!str_detect(family_split, regex("superfamily|subfamily", ignore_case = TRUE))) %>%
    # Remove NAs
    filter(!family_split == "") %>%
    # Count occurrences
    count(family_split, sort = TRUE)
  counts <- counts %>% rename(!!sample := n)
  return(counts)
}

uniprot_gene$Sample <- ifelse(uniprot_gene$Sample == "Dispensable", "Variable", 
                              ifelse(uniprot_gene$Sample == "Private", "Variable", uniprot_gene$Sample))

core_uniprot <- protein_family_function(uniprot_gene[uniprot_gene$Sample == "Core",], "Core")
variable_uniprot <- protein_family_function(uniprot_gene[uniprot_gene$Sample == "Variable",], "Variable")

# merge result to 1 dataframe
merged_families <- core_uniprot %>%
  full_join(variable_uniprot, by = "family_split") %>%
  mutate(across(Core:Variable, ~replace_na(., 0)))
fwrite(merged_families, paste0("uniprot/pan_gene_families.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

# =============================== #
# Statistical Enrichment Analysis #
# =============================== #

# Perform Fisher's Exact Test for Protein Families
# ================================================

# Calculate total core and variable genes across all families
total_core <- sum(merged_families$Core)
total_variable <- sum(merged_families$Variable)

# Initialize vectors to store results
fisher_pvalues <- c()
fisher_odds_ratio <- c()

# Loop over each protein family
for (i in 1:nrow(merged_families)) {
  core_in_family <- merged_families$Core[i]
  variable_in_family <- merged_families$Variable[i]
  core_not_in_family <- total_core - core_in_family
  variable_not_in_family <- total_variable - variable_in_family
  
  # Create contingency table
  contingency_table <- matrix(c(core_in_family, variable_in_family,
                                core_not_in_family, variable_not_in_family),
                              nrow = 2, byrow = TRUE)
  # Fisher's Exact Test
  fisher_test <- fisher.test(contingency_table)
  
  fisher_pvalues <- c(fisher_pvalues, fisher_test$p.value)
  fisher_odds_ratio <- c(fisher_odds_ratio, fisher_test$estimate)
}

# Add results to the dataframe
merged_families <- merged_families %>%
  mutate(
    Fisher_pvalue = fisher_pvalues,
    Fisher_odds_ratio = fisher_odds_ratio
  )

# Adjust p-values for multiple testing using FDR
merged_families <- merged_families %>%
  mutate(Fisher_FDR = p.adjust(Fisher_pvalue, method = "fdr"))
# Interpret enrichment direction and significance
merged_families <- merged_families %>%
  mutate(
    Enrichment = case_when(
      Fisher_odds_ratio > 1 ~ "More Core",
      Fisher_odds_ratio < 1 ~ "More Variable",
      TRUE ~ "No Change"  # in case odds ratio = exactly 1
    ),
    Significance = ifelse(Fisher_FDR < 0.05, "Significant", "Not Significant")
  )
# Save the updated table
fwrite(merged_families, "uniprot/Protein_family_enrichment_results.txt", sep = "\t", quote = F, row.names = F, col.names = T)



# ============= #
# Visualization #
# ============= #
# pangenome categories

plot_protein_fam_significance <- function(df, plot_title,type, output_filename, w,h) {
  
  # get name of the first column
  col <- names(df)[1]
  # Define column order
  column_order <- c("Core", "Variable", "Significance")
  
  # Get top 30 families by Core count
  top_families_by_core <- df %>%
    arrange(desc(Variable)) %>%
    slice_head(n = 30) %>%
    pull(!!sym(col))
  
  # Create significance marks
  df <- df %>%
    mutate(Signif_Mark = case_when(
      Fisher_FDR < 0.01 ~ "**",
      Fisher_FDR  < 0.05 ~ "*",
      TRUE ~ ""
    ))
  
  # Core + Variable long format
  plot_data <- df %>%
    filter(.data[[col]] %in% top_families_by_core) %>%
    mutate(!!col := factor(.data[[col]], levels = rev(top_families_by_core))) %>%
    select(all_of(c(col, "Core", "Variable"))) %>%
    pivot_longer(
      cols = c(Core, Variable),
      names_to = "Category",
      values_to = "Count"
    )  %>%
    mutate(Count = as.character(Count), `Odd's_ratio` = NA_real_)
  
  # Significance column data
  signif_plot_data <- df %>%
    filter(.data[[col]] %in% top_families_by_core) %>%
    mutate(
      !!col := factor(.data[[col]], levels = rev(top_families_by_core)),
      Category = "Significance",
      `Odd's_ratio` = Fisher_odds_ratio,
      Count = Signif_Mark
    ) %>%
    select(col, Category, `Odd's_ratio`, Count)
  
  # Merge
  final_plot_data <- bind_rows(plot_data, signif_plot_data) %>%
    mutate(Category = factor(Category, levels = column_order))
  
  # Plot
  p <- ggplot(final_plot_data, aes(x = Category, y = .data[[col]], fill = `Odd's_ratio`)) +

    geom_tile(color = "black") +
    geom_text(aes(label = Count), size = 3) +
    scale_fill_gradient(
      low = "#92B9C6",
      high = "#B44C32",
      na.value = "white",
      name = "Odds Ratio"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1,angle = 45)) +
    labs(
      x = "Pangenome Category",
      y = type ,
      title = plot_title
    )
  
  # Save as JPEG
  ggsave(filename = paste0(output_filename, ".jpeg"), 
         plot = p, 
         width = w, height = h, 
         units = "in", 
         dpi = 300)
}
plot_protein_fam_significance(merged_families, "Top 30 Protein Families", "Protein Families", "protein_family_enrichment", 7,8) 


# ======= #
# Domains #
# ======= #

# create function for parsing
uniprot_gene_domains_fx <- function(df, sample) {
  counts <- df %>%
  mutate(
    # Extract all domains: anything after "note=" up to the first ";" (non-greedy)
    domain_notes = str_extract_all(`Domain [FT]`, "note=[^;]+")
  ) %>%
  unnest(domain_notes) %>%
  mutate(
    # Remove the "note=" prefix to clean it
    domain_notes = str_remove(domain_notes, "^note="),
    # Remove extra spaces if needed
    domain_notes = str_trim(domain_notes)) %>%
  count(domain_notes, sort = TRUE)
  counts <- counts %>% rename(!!sample := n)
}

core_gene_domains <- uniprot_gene_domains_fx(uniprot_gene[uniprot_gene$Sample == "Core",], "Core")
variable_gene_domains <- uniprot_gene_domains_fx(uniprot_gene[uniprot_gene$Sample == "Variable",], "Variable")
merged_domains <- core_gene_domains %>%
  full_join(variable_gene_domains, by = "domain_notes") %>%
  mutate(across(Core:Variable, ~replace_na(., 0)))


# Calculate total core and variable genes across all families
total_core <- sum(merged_domains$Core)
total_variable <- sum(merged_domains$Variable)

# Initialize vectors to store results
fisher_pvalues <- c()
fisher_odds_ratio <- c()

# Loop over each protein family
for (i in 1:nrow(merged_domains)) {
  core_in_family <- merged_domains$Core[i]
  variable_in_family <- merged_domains$Variable[i]
  core_not_in_family <- total_core - core_in_family
  variable_not_in_family <- total_variable - variable_in_family
  
  # Create contingency table
  contingency_table <- matrix(c(core_in_family, variable_in_family,
                                core_not_in_family, variable_not_in_family),
                              nrow = 2, byrow = TRUE)
  # Fisher's Exact Test
  fisher_test <- fisher.test(contingency_table)
  
  fisher_pvalues <- c(fisher_pvalues, fisher_test$p.value)
  fisher_odds_ratio <- c(fisher_odds_ratio, fisher_test$estimate)
}

# Add results to the dataframe
merged_domains <- merged_domains %>%
  mutate(
    Fisher_pvalue = fisher_pvalues,
    Fisher_odds_ratio = fisher_odds_ratio
  )

# Adjust p-values for multiple testing using FDR
merged_domains <- merged_domains %>%
  mutate(Fisher_FDR = p.adjust(Fisher_pvalue, method = "fdr"))
# Interpret enrichment direction and significance
merged_domains <- merged_domains %>%
  mutate(
    Enrichment = case_when(
      Fisher_odds_ratio > 1 ~ "More Core",
      Fisher_odds_ratio < 1 ~ "More Variable",
      TRUE ~ "No Change"  # in case odds ratio = exactly 1
    ),
    Significance = ifelse(Fisher_FDR < 0.05, "Significant", "Not Significant")
  )
# Save the updated table
fwrite(merged_domains, "uniprot/Protein_domains_enrichment_results.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# ================================================
# visualize using previous function

plot_protein_fam_significance(merged_domains, "Top 30 Protein Domains", "Protein Domains", "protein_domain_enrichment", 5,8) 


# =============== #
# Comparing tests #
# =============== #

install.packages("RVAideMemoire")
library(RVAideMemoire)
# Construct a matrix from merged_families
fam_matrix <- as.matrix(merged_families[, c("Core", "Variable")])
rownames(fam_matrix) <- merged_families$family_split  # or domain_notes for domain test

# Run multi-group Fisher's test with FDR correction
multi_fisher_results <- fisher.multcomp(fam_matrix, p.method = "fdr")
fwrite(multi_fisher_results[["p.value"]], "uniprot/Protein_family_enrichment_results_Raide.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# View the result
print(multi_fisher_results)
