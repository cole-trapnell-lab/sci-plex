# devtools::install_github(repo = 'cole-trapnell-lab/monocle3', ref = "d4a9a35")

# Input path to github directory
path_to_github = "~/sci-Plex/"

# Set directory for sciPlex bin
bin_directory = paste0(path_to_github,
                       "bin",
                       sep = "")

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(viridis)
  library(pheatmap)
  library(scran)
  library(piano)
  library(plyr)
  library(dplyr)
  library(MESS)
  library(tidymodels)
  library(drc)
  library(piano)
  library(snowfall)
  library(monocle3)
  source(paste0(bin_directory,
                "cell_cycle.R",
                sep = ""))
  source(paste0(bin_directory,
                "dose_response.R",
                sep = ""))
  source(paste0(bin_directory,
                "viability.R",
                sep = ""))
  cc.genes = readRDS(paste0(bin_directory,
                            "cc.genes.RDS",
                            sep = ""))
  source(paste0(bin_directory,
                "dispersions_functions.R",
                sep = ""))
  source(paste0(bin_directory,
                "test_drug_dependence.R",
                sep = ""))
})

setwd(paste0(path_to_github,
             "large_screen",
             sep = ""))


# Load previously computed DEG results
A549_diff_test = read.table("Supplementary_Table_5_A549.txt",
                            header = T)


K562_diff_test = read.table("Supplementary_Table_5_K562.txt",
                            header = T)

MCF7_diff_test = read.table("Supplementary_Table_5_MCF7.txt",
                            header = T)


# Pull in drug annotations to annotate DEG results
drug_annotations = read.table("bin/drugProperties_short.tsv",
                              header = T,
                              sep = "\t")

drug_annotations =
  drug_annotations %>%
  dplyr::rename(
    catalog_number = Catalog.Number,
    product_name = Product.Name,
    pathway = Pathway,
    target = Target
  )


# Pull in drug pathway annotations to annotate DEG results
pathway_annotations = read.table("bin/Supplementary_Table_3.txt",
                                 header = T,
                                 sep = "\t")

pathway_annotations =
  pathway_annotations %>%
  dplyr::select(pathway_level_1,
                pathway_level_2,
                catalog_number) %>%
  dplyr::mutate(catalog_number  = as.character(catalog_number)) %>%
  dplyr::distinct()


pathway_annotations$catalog_number[is.na(pathway_annotations$catalog_number)] = "S0000"

drug_annotations = left_join(drug_annotations, pathway_annotations, by = "catalog_number")

diff_test.list = list()
diff_test.list[["A549"]] = left_join(A549_diff_test,
                                     drug_annotations,
                                     by = c("treatment" = "catalog_number"))
diff_test.list[["K562"]] = left_join(K562_diff_test,
                                     drug_annotations,
                                     by = c("treatment" = "catalog_number"))
diff_test.list[["MCF7"]] = left_join(MCF7_diff_test,
                                     drug_annotations,
                                     by = c("treatment" = "catalog_number"))


deg.list = list()


for (cell_line in names(diff_test.list)) {
  diff_test.list[[cell_line]]  = diff_test.list[[cell_line]][grep("dose", diff_test.list[[cell_line]]$term), ]
  
  # Filter genes that have significant dose dependent DEGs
  deg.list[[cell_line]] =
    diff_test.list[[cell_line]] %>%
    filter(q_value < 0.05) %>%
    pull(id) %>%
    unique()
  
  print(length(deg.list[[cell_line]]))
  
}


correlation_matrix.list = list()
correlation_matrix_long.list = list()

for (cell_line in names(diff_test.list)) {
  # Remove genes that have significant viability terms and dose_depdendent terms
  diff_test.list[[cell_line]] = diff_test.list[[cell_line]][(diff_test.list[[cell_line]]$id %in%
                                                               deg.list[[cell_line]]) , ]
  
  diff_test.list[[cell_line]]$id =
    diff_test.list[[cell_line]]$id %>%
    as.character()
  
  diff_test.list[[cell_line]]$treatment =
    diff_test.list[[cell_line]]$treatment %>%
    as.character()
  
  diff_test.list[[cell_line]] =
    diff_test.list[[cell_line]] %>%
    group_by(treatment, id) %>%
    distinct()
  
  # Remove vehicle
  diff_test.list[[cell_line]] = diff_test.list[[cell_line]][diff_test.list[[cell_line]]$product_name != "Vehicle", ]
  
  
  diff_test.list[[cell_line]] =
    dcast(diff_test.list[[cell_line]], id ~ product_name, value.var = "normalized_effect")
  
  row.names(diff_test.list[[cell_line]]) = diff_test.list[[cell_line]]$id
  diff_test.list[[cell_line]]$id = NULL
  
  
  diff_test.list[[cell_line]] = apply(diff_test.list[[cell_line]],
                                      2, as.numeric)
  
  correlation_matrix.list[[cell_line]] = cor(diff_test.list[[cell_line]],
                                             use = "complete.obs", method = "pearson")
  
  
}


# Add pathway annotation for vehicle

ann_colors = list(
  "Pathway" = c(
    "Antioxidant" = "aquamarine3",
    "Apoptotic regulation" = "lemonchiffon4",
    "Cell cycle regulation" = "deepskyblue2",
    "DNA damage & DNA repair" = "slategray4",
    "Epigenetic regulation" = "navy",
    "Focal adhesion signaling" = "brown4",
    "HIF signaling" = "darkgreen",
    "JAK/STAT signaling" = "orangered2",
    "Metabolic regulation" = "gold2",
    "Neuronal signaling" = "darkolivegreen",
    "Nuclear receptor signaling" = "chartreuse3",
    "PKC signaling" = "plum4",
    "Protein folding & Protein degradation" = "darkorchid2",
    "TGF/BMP signaling" = "darkcyan",
    "Tyrosine kinase signaling" = "firebrick1",
    "Other" = "darkgoldenrod4"
  )
)


pheatmap(
  correlation_matrix.list[["A549"]],
  file = "A549_beta_drug_correlation_matrix.png",
  annotation_colors = ann_colors,
  annotation_col = data.frame(
    row.names = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["A549"]]), ]$product_name,
    Pathway = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["A549"]]), ]$pathway_level_1
  ),
  annotation_row = data.frame(
    row.names = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["A549"]]), ]$product_name,
    Pathway = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["A549"]]), ]$pathway_level_1
  ),
  width = 20,
  height = 20,
  color = viridis(25)
)


pheatmap(
  correlation_matrix.list[["K562"]],
  file = "K562_beta_drug_correlation_matrix.png",
  annotation_col = data.frame(
    row.names = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["K562"]]), ]$product_name,
    Pathway = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["K562"]]), ]$pathway_level_1
  ),
  annotation_row = data.frame(
    row.names = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["K562"]]), ]$product_name,
    Pathway = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["K562"]]), ]$pathway_level_1
  ),
  annotation_colors = ann_colors,
  width = 20,
  height = 20,
  color = viridis(25)
)

pheatmap(
  correlation_matrix.list[["MCF7"]],
  file = "MCF7_beta_drug_correlation_matrix.png",
  annotation_col = data.frame(
    row.names = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["MCF7"]]), ]$product_name,
    Pathway = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["MCF7"]]), ]$pathway_level_1
  ),
  annotation_row = data.frame(
    row.names = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["MCF7"]]), ]$product_name,
    Pathway = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["MCF7"]]), ]$pathway_level_1
  ),
  annotation_colors = ann_colors,
  width = 20,
  height = 20,
  color = viridis(25)
)

correlation_matrix_bounded.list = correlation_matrix.list

for (cell_line in names(correlation_matrix_bounded.list)) {
  correlation_matrix_bounded.list[[cell_line]][correlation_matrix_bounded.list[[cell_line]] > 0.6] = 0.6
  
}

pheatmap(
  correlation_matrix_bounded.list[["A549"]],
  file = "supplemental_figure15_A549_beta_drug_correlation_matrix_bounded.png",
  annotation_col = data.frame(
    row.names = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["A549"]]), ]$product_name,
    Pathway = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["A549"]]), ]$pathway_level_1
  ),
  annotation_row = data.frame(
    row.names = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["A549"]]), ]$product_name,
    Pathway = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["A549"]]), ]$pathway_level_1
  ),
  annotation_colors = ann_colors,
  width = 10,
  height = 10,
  fontsize = 3.5,
  annotation_legend = F,
  legend = F,
  treeheight_row = 3,
  treeheight_col = 3,
  color = viridis(25)
)


pheatmap(
  correlation_matrix_bounded.list[["K562"]],
  file = "supplemental_figure16_K562_beta_drug_correlation_matrix_bounded.png",
  annotation_col = data.frame(
    row.names = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["K562"]]), ]$product_name,
    Pathway = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["K562"]]), ]$pathway_level_1
  ),
  annotation_row = data.frame(
    row.names = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["K562"]]), ]$product_name,
    Pathway = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["K562"]]), ]$pathway_level_1
  ),
  annotation_colors = ann_colors,
  width = 10,
  height = 10,
  fontsize = 3.5,
  annotation_legend = F,
  legend = F,
  treeheight_row = 3,
  treeheight_col = 3,
  color = viridis(25)
)

pheatmap(
  correlation_matrix_bounded.list[["MCF7"]],
  file = "supplemental_figure17_MCF7_beta_drug_correlation_matrix_bounded.png",
  annotation_col = data.frame(
    row.names = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["MCF7"]]), ]$product_name,
    Pathway = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["MCF7"]]), ]$pathway_level_1
  ),
  annotation_row = data.frame(
    row.names = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["MCF7"]]), ]$product_name,
    Pathway = drug_annotations[drug_annotations$product_name %in% colnames(correlation_matrix.list[["MCF7"]]), ]$pathway_level_1
  ),
  annotation_colors = ann_colors,
  width = 10,
  height = 10,
  fontsize = 3.5,
  annotation_legend = F,
  legend = F,
  treeheight_row = 3,
  treeheight_col = 3,
  color = viridis(25)
)