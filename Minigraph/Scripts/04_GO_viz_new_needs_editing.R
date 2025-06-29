source("/Users/jongpaduhilao/Desktop/LAB_Files/Initial_Pangenome_analysis/Trial_4/Scripts/04_GO_viz.R")

plot_GO <- function(go_mutated,title) {
  go_mutated <- read.csv(go_mutated, header = TRUE)
  
  accession_order <- unique(go_mutated$accession)
  go_mutated$accession <- factor(go_mutated$accession, levels = accession_order)
  
  plt <- ggplot(go_mutated, aes(x = accession, y = `GO_terms`)) +
    geom_point(aes(size = `gene_count`, color = accession)) +
    scale_size_continuous(range = c(1, 10)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    labs(x = "Accession", y = "Gene Ontology Terms", size = "Number of genes",title = title ) +
    guides(color = FALSE)
  
  print(plt)
}

plot_GO_over <- function(go_mutated,title, top_n = 10) {
  go_mutated <- read.csv(go_mutated, header = TRUE)
  # Ensure FDR is numeric and log-transform it
  go_mutated$FDR <- as.numeric(go_mutated$FDR)
  go_mutated$log_FDR <- -log10(go_mutated$FDR)
  
  first_column_name <- names(go_mutated)[1]
  accession_order <- unique(go_mutated$Accession)
  
  # filter tops
  go_filtered <- go_mutated %>%
    group_by(Accession) %>%
    slice_min(order_by = FDR, n = top_n) %>%
    ungroup()
  
  go_filtered$Accession <- factor(go_filtered$Accession, levels = accession_order)
  
  # plot
  
  plt <- ggplot(go_filtered, aes(x = Accession, y = .data[[first_column_name]])) + # edit Y 
    geom_point(aes(size = Query, color = log_FDR)) +  
    # scale_shape_manual(values = c("+" = 16, "-" = 17)) +
    scale_size_continuous(range = c(1, 9)) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(FDR)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    labs(x = "Accession", y = "Gene Ontology Terms", 
         size = "Number of genes", 
         title = title) +
    guides(color = guide_colorbar(barwidth = 1, barheight = 3))
    
    print(plt)
}

GO_mf <- "/Users/jongpaduhilao/Downloads/Non_reference_GO/unfiltered_allOrg_MF.csv"
GO_pc <- "/Users/jongpaduhilao/Downloads/Non_reference_GO/unfiltered_allOrg_PC.csv"

bp_fdr <- "/Users/jongpaduhilao/Desktop/LAB_Files/Initial_Pangenome_analysis/Trial_4/03_Non_reference_GO/unfiltered_enrichment_BP.csv"
pc_fdr <- "/Users/jongpaduhilao/Desktop/LAB_Files/Initial_Pangenome_analysis/Trial_4/03_Non_reference_GO/unfiltered_enrichment_PC.csv"
mf_fdr <- "/Users/jongpaduhilao/Desktop/LAB_Files/Initial_Pangenome_analysis/Trial_4/03_Non_reference_GO/unfiltered_enrichment_MF.csv"


pc_fdr <- "/Users/jongpaduhilao/Downloads/Non_reference_GO/unfiltered_enrichment_PC.csv"
mf_fdr <- "/Users/jongpaduhilao/Downloads/Non_reference_GO/unfiltered_enrichment_MF.csv"

mf <- plot_GO(GO_mf,"Molecular Function(MF)")
pc <- plot_GO(GO_pc, "Protein Class(PC)") 


##### GO with P-values #####
cw02_fdr <- "/Users/jongpaduhilao/Downloads/Non_reference_GO/unfiltered_enrichment_PC.csv"
go_mutated <- read.csv(cw02_fdr, header = TRUE)
pc_enrich <- plot_GO_over(pc_fdr, "Protein Class Enrichment",10)
mf_enrich <- plot_GO_over(mf_fdr, "Molecular Function Enrichment",10)
bp_enrich <- plot_GO_over(bp_fdr, "Biological Process Enrichment", 10)

##### GO for reference blast hits #####

ref_bp <- "/Users/jongpaduhilao/Desktop/LAB_Files/Initial_Pangenome_analysis/Trial_4/03_Core_GO/BP_enrich_categories.csv"
ref_mf <-"/Users/jongpaduhilao/Desktop/LAB_Files/Initial_Pangenome_analysis/Trial_4/03_Core_GO/MF_enrich_categories.csv"
ref_pc <- "/Users/jongpaduhilao/Desktop/LAB_Files/Initial_Pangenome_analysis/Trial_4/03_Core_GO/PC_enrich_categories.csv"

ref_bp_plot <- plot_GO_over(ref_bp, "Biological Process Enrichment", 10)
ref_mf_plot <- plot_GO_over(ref_mf, "Molecular Function Enrichment", 10)
ref_pc_plot <- plot_GO_over(ref_pc, "Protein Class Enrichment", 10)


