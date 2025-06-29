library(ggplot2)
library(reshape2)
library(dplyr)
source("../../Scripts/03_Ref_Analysis.R")
source("../../Scripts/04_NonRef_Analysis.R")

PA_matrix <- fread("02_Presence_Absence_ALL.tsv")

# to categorize the pangenome on-reference, use script 03_Ref_Analysis.R
reference_categorization(PA_matrix)

# for plotting, use script 03_Ref_Analysis.R
# plot has some boxes, need to reset ggplot theme here
intersect_ref <- intersection(PA_matrix)

theme_set(theme_bw(base_size = 9)+
            theme(panel.grid = element_blank()))
upset(fromList(intersect_ref), nsets=6, nintersects = 40, order.by = c("freq"), 
      set_size.show = FALSE,  mainbar.y.label = "Segment Number", sets.x.label = "Segment Number")

# Use script 03_Ref_Analysis.R to quantify the length of sharing pattern length
plt_ref <- reference_sharing(PA_matrix)
print(plt_ref)

#####################terminal###############################
## In server, acquire the core, dispensable, and private tsv files
## Grab IDs using the terminal and output as txt for downstream analysis
## Use the terminal commands for more information 
###########################terminal##########################

############# GO visualization ##############################
# After getting IDs and Gene Ontology enrichment, we can now visualize the data here. 
colsNames <-  c("number","Gene Ontology", "Number of genes",
                "%gene hits against total gene hits",
                "%gene hits against total function hits")
## Biological Process
bp_core <- fread("03_Gene_Ontology/functional_class_BP_core.txt")
colnames(bp_core) <- colsNames
bp_core$category <- "core"

bp_disp <- fread("03_Gene_Ontology/functional_class_BP_disp.txt")
colnames(bp_disp) <- colsNames
bp_disp$category <- "dispensable"

bp_priv <- fread("03_Gene_Ontology/functional_class_BP_private.txt")
colnames(bp_priv) <- colsNames
bp_priv$category <- "private"

bp_all <- rbind(bp_core, bp_disp, bp_priv)

## Molecular function
mf_core <- fread("03_Gene_Ontology/functional_class_MF_core.txt")
colnames(mf_core) <- colsNames
mf_core$category <- "core"

mf_disp <- fread("03_Gene_Ontology/functional_class_MF_disp.txt")
colnames(mf_disp) <- colsNames
mf_disp$category <- "dispensable"

mf_private <- fread("03_Gene_Ontology/functional_class_MF_private.txt")
colnames(mf_private) <- colsNames
mf_private$category <- "private"

mf_all <- rbind(mf_core, mf_disp, mf_private)

## Protein class
prot_core <- fread("03_Gene_Ontology/functional_class_protein_core.txt")
colnames(prot_core) <- colsNames
prot_core$category <- "core"

prot_disp <- fread("03_Gene_Ontology/functional_class_protein_disp.txt")
colnames(prot_disp) <- colsNames
prot_disp$category <- "dispensable"

prot_private <- fread("03_Gene_Ontology/functional_class_protein_private.txt")
colnames(prot_private) <- colsNames
prot_private$category <- "private"

prot_all <- rbind(prot_core, prot_disp, prot_private)

### visualize
source("/Users/jongpaduhilao/Desktop/LAB Files/Initial_Pangenome_analysis/Trial_4/Scripts/04_GO_viz.R")
bp_viz <- go_visualize_1(bp_all)
print(bp_viz)

mf_viz <- go_visualize_1(mf_all)
print(mf_viz)

prot_viz <- go_visualize_1(prot_all)
print(prot_viz)

###### ALL Privates and combinations###############
prot_IRGSP <- fread("03_Gene_Ontology_private/functional_class_protein_IRGSP.txt")
colnames(prot_IRGSP) <- colsNames
prot_IRGSP$category <- "IRGSP"

prot_nh232 <- fread("03_Gene_Ontology_private/functional_class_protein_nh232.txt")
colnames(prot_nh232) <- colsNames
prot_nh232$category <- "nh232"

prot_cw02 <- fread("03_Gene_Ontology_private/functional_class_protein_cw02.txt")
colnames(prot_cw02) <- colsNames
prot_cw02$category <- "cw02"

prot_nh236 <- fread("03_Gene_Ontology_private/functional_class_protein_nh236.txt")
colnames(prot_nh236) <- colsNames
prot_nh236$category <- "nh236"

prot_nh286 <- fread("03_Gene_Ontology_private/functional_class_protein_nh286.txt")
colnames(prot_nh286) <- colsNames
prot_nh286$category <- "nh286"

prot_nh273 <- fread("03_Gene_Ontology_private/functional_class_protein_nh273.txt")
colnames(prot_nh273) <- colsNames
prot_nh273$category <- "nh273"

prot_asian <- fread("03_Gene_Ontology/functional_class_protein_asian.txt")
colnames(prot_asian) <- colsNames
prot_asian$category <- "asian rice only"

prot_african <- fread("03_Gene_Ontology/functional_class_protein_african.txt")
colnames(prot_african) <- colsNames
prot_african$category <- "african rice only"

prot_domesticated <- fread("03_Gene_Ontology/functional_class_protein_domesticated.txt")
colnames(prot_domesticated) <- colsNames
prot_domesticated$category <- "domesticated only"

prot_wild <- fread("03_Gene_Ontology/functional_class_protein_wild.txt")
colnames(prot_wild) <- colsNames
prot_wild$category <- "wild species only"

### combine all to one dataframe
prot_priv_all <- rbind(prot_asian, prot_african, prot_domesticated, prot_wild,
                       prot_IRGSP, prot_nh232, prot_cw02, prot_nh236, prot_nh286, prot_nh273)

### visualize
prot_priv_viz <- go_visualize_2(prot_priv_all)
print(prot_priv_viz)
