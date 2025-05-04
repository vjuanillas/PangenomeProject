library(data.table)
library(dplyr)
library(tidyverse)
library(car)            # For leveneTest
library(rstatix)        # For kruskal_test and post-hoc tests
library(ggpubr)         # For nice visualizations

# set wd
# change the path to vcf if needed
setwd("/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p60/")
# ==== #
# data #
# ==== #
vcf_clean <- fread("whole_genome/VCF_Parsimony_Not_raw_statistics.txt")
vcf_parsimonious <- fread("whole_genome/VCF_Clean_parsimony_labelled_gain_loss_final.txt")
# ================ #
# Statistical Test #
# ================ #
# 1. Non-parsimonious cases
# 2. parsimonious cases

if (!dir.exists("statistical_tests")) {
  dir.create("statistical_tests")
}

# ====================== #
# Non-parsimonious cases #
# ====================== #
# Question: Is there a significant difference between the different non-parsimonious cases in terms of:
# a. Frequency and b. Cumulative Length


# View Summary of Non-parsimonious cases
vcf_non_parsimonious <- vcf_clean %>% filter(par_not == "not")
vcf_non_parsimonious_stats <- vcf_non_parsimonious[, .(
  N_bubble = .N, 
  REF_length = sum(nchar(REF)),
  ALT_length = sum(nchar(ALT))
), by = .(ALT_accn, type, subtype)] %>% 
  mutate(`Variant Length` = ifelse(subtype %in% c("INS", "SV_INS"), ALT_length, REF_length)) %>%
  filter(subtype %in% c("INS","DEL","SV_INS","SV_DEL"))

vcf_non_parsimonious_stats_grouped <- vcf_non_parsimonious_stats[, .(
  Count = sum(N_bubble),
  Variant_length = sum(`Variant Length`),
  Average_length = sum(`Variant Length`)/sum(N_bubble)
), by =.(ALT_accn, type)]


### Conduct Test: COUNT
# Separate by type
types <- unique(vcf_non_parsimonious_stats_grouped$type)

for (t in types) {
  cat("\n====== Testing for type:", t, "======\n")
  
  df_type <- vcf_non_parsimonious_stats_grouped %>% filter(type == t)
  
  ## 1. Frequency/Count testing (Chi-squared)
  cat("\n- Chi-square test on Counts:\n")
  chisq_counts <- chisq.test(df_type$Count)
  print(chisq_counts)
  residuals_df <- data.frame(
    ALT_accn = df_type$ALT_accn,
    residuals = chisq_counts$stdres
  )
  fwrite(residuals_df, paste0("statistical_tests/Chisq_count_residuals_",t,".txt"), sep = "\t", 
         quote=F, row.names = F, col.names=T)
  
  plot <- ggplot(residuals_df, aes(x = ALT_accn, y = residuals, fill = residuals > 0)) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("#B44C32", "#92B9C6"), guide = "none") +
    labs(
      title = paste0("Standardized Residuals from Chi-square Test","(",t,")"),
      x = "Non-parsimonious Cases",
      y = "Standardized Residuals"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  ggsave(
    filename = paste0("statistical_tests/Chisq_count_residuals_",t,".jpeg"),
    plot = plot,
    width = 11,
    height = 9,
    device = "jpeg",
    dpi = 300
  )
}

### Conduct Test: Variant Length
for (t in types) {
  cat("\n====== Testing for type:", t, "======\n")
  df_type <- vcf_non_parsimonious_stats_grouped %>% filter(type == t)  
  ## 2. Variant Length testing
  cat("\n- Checking normality of Variant Lengths (Shapiro-Wilk):\n")
  shapiro <- shapiro.test(df_type$Variant_length)
  print(shapiro)
  
  cat("\n- Checking equal variance (Levene’s Test):\n")
  levene <- leveneTest(Variant_length ~ ALT_accn, data = df_type)
  print(levene)
  
  ## 3. Decide which test: ANOVA or Kruskal
  if (shapiro$p.value > 0.05 & levene$`Pr(>F)`[1] > 0.05) {
    cat("\n- Data is normal and variances are equal: Using ANOVA\n")
    aov_res <- aov(Variant_length ~ ALT_accn, data = df_type)
    summary(aov_res)
    cat("\n- Post-hoc Tukey test:\n")
    print(TukeyHSD(aov_res))
  } else {
    cat("\n- Data is NOT normal or variances unequal: Using Kruskal-Wallis\n")
    kruskal <- kruskal_test(Variant_length ~ ALT_accn, data = df_type)
    print(kruskal)
    
    cat("\n- Post-hoc Dunn's test:\n")
    dunn_res <- df_type %>%
      dunn_test(Variant_length ~ ALT_accn)
    print(dunn_res)
  }
}

# ================== #
# Parsimonious cases #
# ================== #
library(purrr)
library(broom) # for clean outputs
# Question: Is there a significant difference between the different parsimonious cases within INDELs and SVs in terms of:
# a. Frequency and b. Cumulative Length
vcf_parsimonious_INDEL <- vcf_parsimonious %>% filter(type == "INDEL")
vcf_parsimonious_SV <- vcf_parsimonious %>% filter(type == "SV")

vcf_parsimonious_INDEL_stats <- vcf_parsimonious_INDEL[, .(
  N_bubble = .N,
  Variant_length = sum(N_bp_diff)
), by = .(Taxon, Gain_Loss)]

vcf_parsimonious_SV_stats <- vcf_parsimonious_SV[, .(
  N_bubble = .N,
  Variant_length = sum(N_bp_diff)
), by = .(Taxon, Gain_Loss)]

# FOR INDEL counts
# Perform chi-square tests for each taxon in a single pipeline
chi_results_fx <- function(df) {
  chi_results <- df  %>%
  # Group by taxon
  group_by(Taxon) %>%
  # Perform chi-square test for each taxon
  summarize(
    chi_test = list(chisq.test(N_bubble, p = rep(0.5, n()))),
    .groups = "drop"
  ) %>%
  # Extract and tidy results
  mutate(test_results = map(chi_test, tidy)) %>%
  unnest(test_results) %>%
  # Clean up output
  select(Taxon, statistic, p.value, parameter, method)
  
  # Add FDR correction
  chi_results_with_fdr <- chi_results %>%
  mutate(
    # Add FDR correction (Benjamini-Hochberg method)
    p.value.adjusted = p.adjust(p.value, method = "fdr"),
    # Format p-values for better readability
    p.value = format.pval(p.value, digits = 3),
    p.value.adjusted = format.pval(p.value.adjusted, digits = 3),
    significant = p.value.adjusted < 0.05
  )
  print(chi_results_with_fdr)
}
chi_results_fx(vcf_parsimonious_INDEL_stats)
chi_results_fx(vcf_parsimonious_SV_stats)$residuals
