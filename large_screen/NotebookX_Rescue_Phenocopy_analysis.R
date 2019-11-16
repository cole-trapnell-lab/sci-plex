# devtools::install_github(repo = 'cole-trapnell-lab/monocle3', ref = "d4a9a35")

# Input path to github directory
path_to_github = "~/sci-Plex/"

# Set directory for sciPlex bin
bin_directory = paste0(path_to_github,
                       "bin",
                       sep = "")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(tidyr)
  library(Matrix)
  library(scran)
  library(pheatmap)
  library(monocle3)
})

# Set DelayedArray Parameters
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)


# Load sci-Chem specific scripts
suppressPackageStartupMessages({
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

###############
# FIXME : MNN function call put in the script that calls that 
###############

# User functions
refresh_cds = function(cds){
  mat = counts(cds)
  coldata =
    colData(cds) %>%
    as.data.frame()
  rowdata =
    rowData(cds) %>%
    as.data.frame()

  cds = new_cell_data_set(expression_data = mat,
                          cell_metadata = coldata,
                          gene_metadata = rowdata)
}

setwd(paste0(path_to_github,
             "large_screen",
             sep = ""))

#### Load cds object of full 188 compound screen and subset by human genes and cells treated with HDAC inhibitors
screen_cds <- readRDS("cds.RDS")
screen_hdac_cells <- read.table("hdac_cells.txt")
screen_hdac_cells <- as.character(screen_hdac_cells$V1)

# Subset human genes
human.genes = rowData(screen_cds)$id[grepl("^ENSG", rowData(screen_cds)$id)]

screen_hdac_cds <- screen_cds[rowData(screen_cds)$id %in% human.genes,
                              colData(screen_cds)$cell %in% screen_hdac_cells]

rm(screen_cds)

screen_hdac_cds <- screen_hdac_cds[,!is.na(colData(screen_hdac_cds)$time_point)]

# Re-work columns to match rescue phenocopy experiment
pData(screen_hdac_cds) <- pData(screen_hdac_cds)[,c("cell", "sample","Size_Factor","n.umi",
                                      "cell_type","product_name","dose")]

colnames(pData(screen_hdac_cds)) <- c("cell", "sample","Size_Factor","n.umi",
                                      "cell_type","hdac","hdac_dose")

pData(screen_hdac_cds)$cell_type <- as.character(pData(screen_hdac_cds)$cell_type)
pData(screen_hdac_cds)$hdac <- as.character(pData(screen_hdac_cds)$hdac)

#### Load cds for rescue/phenocopy experiment
cds_hdaci_mod = readRDS("sciPlex_4.RDS")

cell_meta_data = stringr::str_split_fixed(colData(cds_hdaci_mod)$Cell,pattern = "_",n = 8)

colData(cds_hdaci_mod)$pcr_well = cell_meta_data[,1] %>% as.character()
colData(cds_hdaci_mod)$pcr_plate = cell_meta_data[,2] %>% as.character()
colData(cds_hdaci_mod)$rt_well = cell_meta_data[,5] %>% as.character()
colData(cds_hdaci_mod)$lig_well = cell_meta_data[,8] %>% as.character()

hash_meta_data = stringr::str_split_fixed(colData(cds_hdaci_mod)$top_oligo,pattern = "_",n = 7)

colData(cds_hdaci_mod)$treatment_dose = hash_meta_data[,1] %>% as.numeric()
colData(cds_hdaci_mod)$treatment = hash_meta_data[,2] %>% as.character()
colData(cds_hdaci_mod)$hdac_dose = hash_meta_data[,3] %>% as.numeric()
colData(cds_hdaci_mod)$hdac = hash_meta_data[,4] %>% as.character()
colData(cds_hdaci_mod)$cell_type = hash_meta_data[,5] %>% as.character()
colData(cds_hdaci_mod)$hash_plate = hash_meta_data[,6] %>% as.character()
colData(cds_hdaci_mod)$hash_well = hash_meta_data[,7] %>% as.character()

validly_labeled_cells =
  colData(cds_hdaci_mod) %>%
  as.data.frame() %>%
  filter(!is.na(top_oligo)) %>%
  filter(hash_umis >= 5,
         top_to_second_best_ratio >= 5) %>%
  pull(Cell)

cds_hdaci_mod = cds_hdaci_mod[,validly_labeled_cells]

pData(cds_hdaci_mod) <- pData(cds_hdaci_mod)[,c("Cell","sample","Size_Factor","n.umi",
                                                "cell_type","hdac","hdac_dose",
                                                "treatment","treatment_dose")]

colnames(pData(cds_hdaci_mod)) <- c("cell", "sample","Size_Factor","n.umi",
                                    "cell_type","hdac","hdac_dose","treatment","treatment_dose")

pData(cds_hdaci_mod)$hdac_dose <- pData(cds_hdaci_mod)$hdac_dose*1000

# Add place holder columns to HDAC inhibitor screen cds to match rescue/phenocopy experiment
pData(screen_hdac_cds)$treatment <- rep(NA, length(row.names(pData(screen_hdac_cds))))
pData(screen_hdac_cds)$treatment_dose <- rep(NA, length(row.names(pData(screen_hdac_cds))))

identical(colnames(pData(cds_hdaci_mod)), colnames(pData(screen_hdac_cds)))
identical(row.names(exprs(cds_hdaci_mod)),row.names(exprs(screen_hdac_cds)))
identical(row.names(rowData(cds_hdaci_mod)), row.names(rowData(screen_hdac_cds)))

# Create a joint cds, useful for cleaning up the experiment for misassigned hashes, for aligned comparisons and is the starting cds for DEG testing
joint_exprs <- cbind(exprs(screen_hdac_cds),exprs(cds_hdaci_mod))
joint_pData <- rbind(colData(screen_hdac_cds),colData(cds_hdaci_mod))
joint_fData <- rowData(screen_hdac_cds)

rm(screen_hdac_cds,cds_hdaci_mod)

identical(colnames(joint_exprs), row.names(joint_pData))
identical(row.names(joint_exprs), row.names(joint_fData))

joint_cds <- new_cell_data_set(expression_data = joint_exprs,
                               cell_metadata = joint_pData,
                               gene_metadata = joint_fData)

rm(joint_exprs, joint_pData, joint_fData)


new_hdac <- sapply(pData(joint_cds)$hdac, function(x){
  
  new_hdac <- strsplit(x, "\\s+")[[1]][1]
  return(new_hdac)
  
})

new_hdac <- sapply(new_hdac, function(x){
  
  new_hdac <- ifelse(x == "DMSO","Vehicle",x)
  return(new_hdac)

})

pData(joint_cds)$new_hdac <- new_hdac

# Process joint cds and perform dimensionality reduction, filter misassigned hash cells
joint_cds <- joint_cds[,pData(joint_cds)$cell_type %in% c("A549","MCF7")]

joint_cds = detect_genes(joint_cds)
joint_cds = estimate_size_factors(joint_cds)

joint_cds <- estimate_cell_cycle(joint_cds,
                                 g1s_markers = cc.genes$s.genes,
                                 g2m_markers = cc.genes$g2m.genes)

genes_expressed_per_drug = colData(joint_cds) %>%
                                as.data.frame() %>%
                                group_by(treatment, hdac, cell_type) %>%
                                nest() %>%
                                mutate(fraction_genes = purrr::map(data, .f = function(pdata_subset, cds) {
                                    cds_subset = cds[,as.character(pdata_subset$cell)]
                                    cds_subset = detect_genes(cds_subset)
                                    tibble(id = rowData(cds_subset)$id,
                                           gene_short_name = rowData(cds_subset)$gene_short_name,
                                           num_cells_expressed = rowData(cds_subset)$num_cells_expressed,
                                           fraction_cells_expressed = rowData(cds_subset)$num_cells_expressed / ncol(cds_subset))
                                }, joint_cds))

genes_expressed_per_drug = genes_expressed_per_drug %>% 
  unnest(fraction_genes)

expressed_genes = genes_expressed_per_drug %>% 
  filter(fraction_cells_expressed > 0.01) %>% 
  ungroup() %>% 
  dplyr::select(id) %>% 
  distinct()

expressed_genes = as.character(expressed_genes$id)
length(expressed_genes)

joint_cds = preprocess_cds(joint_cds,
                           method = "PCA",
                           num_dim = 35,
                           norm_method = "log",
                           residual_model_formula_str = "~log(n.umi)",
                           use_genes = expressed_genes)
                               
joint_cds = reduce_dimension(cds = joint_cds,
                             max_components = 2,
                             reduction_method = "UMAP",
                             umap.metric = "cosine",
                             umap.n_neighbors = 50,
                             umap.min_dist = 0.5,
                             umap.fast_sgd=FALSE, 
                             cores=1,
                             verbose = T)

colData(joint_cds)$umap1 <- joint_cds@reducedDims[["UMAP"]][,1]
colData(joint_cds)$umap2 <- joint_cds@reducedDims[["UMAP"]][,2]

joint_cds = cluster_cells(cds = joint_cds,
                              reduction_method = "UMAP")
                                          
colData(joint_cds)$partitions = joint_cds@clusters[["UMAP"]]$partitions

ggplot(as.data.frame(colData(joint_cds)), aes(x = umap1, y = umap2, color = partitions)) +
geom_point() +
facet_wrap(~cell_type, ncol = 2) +
monocle3:::monocle_theme_opts() + ggsave("test_3.png")#+
#ggsave("figures/UMAP_by_cell_type.png")     

# Remove  the small number of cells that have been missasigned by the hash
# New filters
joint_cds <- joint_cds[,(pData(joint_cds)$partitions == "2" & pData(joint_cds)$cell_type == "A549") |
                          (pData(joint_cds)$partitions == "1" & pData(joint_cds)$cell_type == "MCF7")]

joint_cds <- refresh_cds(joint_cds)

# saveRDS(joint_cds, "Joint_HDAC_cds_clean.RDS")
##### This cds is loaded for DEG testing in a separate notebook #####
# joint_cds <- readRDS("Joint_HDAC_cds_clean.RDS")

# Isolate cells from rescue experiment to independently determine the effect of HDACi's on HDACi pseudotime progression 
cds.list <- list()

cds.list[["A549"]] <- joint_cds[,pData(joint_cds)$cell_type == "A549" &
                                  pData(joint_cds)$sample == "sciChem20"]
cds.list[["MCF7"]] <- joint_cds[,pData(joint_cds)$cell_type == "MCF7" &
                                  pData(joint_cds)$sample == "sciChem20"]

dispersion_test.list <- list()
disp_table.list <- list()

for(cell_line in names(cds.list)){

  dispersion_test.list[[cell_line]] <- estimateDispersionsForCellDataSet(cds.list[[cell_line]], 
                                                                         min_cells_detected = 100, 
                                                                         removeOutliers = FALSE)
  disp_table.list[[cell_line]] <- dispersionTable(dispersion_test.list[[cell_line]])

}

unsup_clustering_genes.list <- list()

for(cell_line in names(cds.list)){
  
  disp_table.list[[cell_line]] <- disp_table.list[[cell_line]] %>% 
  mutate(excess_disp = (dispersion_empirical - dispersion_fit) / dispersion_fit) %>% 
  dplyr::arrange(plyr::desc(excess_disp))
  
  unsup_clustering_genes.list[[cell_line]] <- as.character(head(disp_table.list[[cell_line]], 1000)$gene_id)

}

for(cell_line in names(cds.list)){

cds.list[[cell_line]] = detect_genes(cds.list[[cell_line]])
cds.list[[cell_line]] = estimate_size_factors(cds.list[[cell_line]])

cds.list[[cell_line]] <- estimate_cell_cycle(cds.list[[cell_line]],
                                             g1s_markers = cc.genes$s.genes,
                                             g2m_markers = cc.genes$g2m.genes)

cds.list[[cell_line]] = preprocess_cds(cds.list[[cell_line]],
                                           method = "PCA",
                                           num_dim = 35,
                                           norm_method = "log",
                                           use_genes = unsup_clustering_genes.list[[cell_line]])
  
cds.list[[cell_line]] = reduce_dimension(cds = cds.list[[cell_line]],
                             max_components = 2,
                             reduction_method = "UMAP",
                             umap.metric = "cosine",
                             umap.n_neighbors = 50,
                             umap.min_dist = 0.5,
                             umap.fast_sgd=FALSE, 
                             cores=1,
                             verbose = T)

cds.list[[cell_line]] = cluster_cells(cds = cds.list[[cell_line]],
                          reduction_method = "UMAP",
                          resolution = 1e-6)

colData(cds.list[[cell_line]])$umap1 <- cds.list[[cell_line]]@reducedDims[["UMAP"]][,1]
colData(cds.list[[cell_line]])$umap2 <- cds.list[[cell_line]]@reducedDims[["UMAP"]][,2]

}


for(cell_line in names(cds.list)){

graph_parameters = list()
graph_parameters[["minimal_branch_len"]] = 10
graph_parameters[["ncenter"]] = 100

cds.list[[cell_line]] <- learn_graph(cds.list[[cell_line]], 
                                     learn_graph_control = graph_parameters, 
                                     use_partition = F, 
                                     close_loop = FALSE)

}

for(cell_line in names(cds.list)){
  
  pData(cds.list[[cell_line]])$closest_vertex <- cds.list[[cell_line]]@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[,1]
  
  ordering_summary = 
    pData(cds.list[[cell_line]]) %>% 
    as.data.frame() %>%
    dplyr::group_by(closest_vertex) %>% 
    dplyr::count(hdac_dose) 
  
  ordering_summary = 
    ordering_summary %>% 
    mutate(total_cells = sum(n), fraction_dose = n / total_cells)
  
  root_nodes =  ordering_summary %>% filter(hdac_dose == 0 & fraction_dose > 0.5)
  root_nodes = root_nodes$closest_vertex

  pData(cds.list[[cell_line]])$root_node = pData(cds.list[[cell_line]])$closest_vertex %in% root_nodes
  
  root_cells  = 
    colData(cds.list[[cell_line]]) %>% 
    as.data.frame() %>% 
    filter(root_node) %>% pull(cell) %>% 
    as.character()
  
  cds.list[[cell_line]] <- order_cells(cds.list[[cell_line]], root_cells=root_cells)
  
  pData(cds.list[[cell_line]])$Pseudotime = cds.list[[cell_line]]@principal_graph_aux[["UMAP"]]$pseudotime
  pData(cds.list[[cell_line]])$Pseudotime_bin = cut(pData(cds.list[[cell_line]])$Pseudotime, 
                                                    breaks=quantile(pData(cds.list[[cell_line]])$Pseudotime, 
                                                                    seq(0, 1, 0.1)), 
                                                    labels=F)
  
}

for(cell_line in names(cds.list)){
  
  pData(cds.list[[cell_line]])$dose_character <- as.character(pData(cds.list[[cell_line]])$hdac_dose)
  pData(cds.list[[cell_line]])$dose_character <- factor(pData(cds.list[[cell_line]])$dose_character, 
                                                        levels = c("0", "10", "50", "100", "500",
                                                                   "1000", "5000", "10000"))
}

for(cell_line in names(cds.list)){
  
  pData(cds.list[[cell_line]])$treatment_dose_character <- as.character(pData(cds.list[[cell_line]])$treatment_dose)
  pData(cds.list[[cell_line]])$treatment_dose_character[is.na(pData(cds.list[[cell_line]])$treatment_dose_character)] <- "untreated"
  pData(cds.list[[cell_line]])$treatment_dose_character <- factor(pData(cds.list[[cell_line]])$treatment_dose_character,
                                                                  levels = c("0","0.05","0.25","0.5","1","10","100",
                                                                             "2","2.5","20","25","5","50", "untreated"))
  }

# saveRDS(cds.list, "A549_MCF7_rescue_phenocopy_processed_cds.rds")
# cds.list <- readRDS("A549_MCF7_rescue_phenocopy_processed_cds.rds")

plot_cells(cds.list[["A549"]], reduction_method = "UMAP", color_cells_by = "dose_character", label_cell_groups = FALSE,
           label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE, cell_size = 0.5, trajectory_graph_segment_size = 0.5) +
  scale_color_manual("Dose",values = c("0"="gray", "10"="#1D1147FF", "50"="#51127CFF", "100"="#822681FF", "500"="#B63679FF",
                                       "1000"="#E65164FF", "5000" = "#FB8861FF", "10000"="#FEC287FF")) +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.25,"line")) +
  guides(guides(colour = guide_legend(override.aes = list(size=1)))) + xlim(-10,10) + ylim(-6,6) +
  ggsave("supplemental_figure31_Aligned_A549_HDACi_UMAP_screen_sciplex_pheno_rescue.png", 
         dpi = 600, height = 1.5, width = 1.8)

plot_cells(cds.list[["MCF7"]], reduction_method = "UMAP", color_cells_by = "dose_character", label_cell_groups = FALSE,
           label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE, cell_size = 0.5, trajectory_graph_segment_size = 0.5) +
  #facet_wrap(~treatment) +
  scale_color_manual("Dose",values = c("0"="gray", "10"="#1D1147FF", "50"="#51127CFF", "100"="#822681FF", "500"="#B63679FF",
                                       "1000"="#E65164FF", "5000" = "#FB8861FF", "10000"="#FEC287FF")) +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.25,"line")) #+
  guides(guides(colour = guide_legend(override.aes = list(size=1)))) +
  ggsave("supplemental_figure31_Aligned_MCF7_HDACi_UMAP_screen_sciplex_pheno_rescue.png", 
         dpi = 600, height = 1.5, width = 1.8)


colors <- colorRampPalette(c("midnightblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1"))(10)

# Make a copy of the cds for plotting where bins are in order
# TO DO: change to just copy pData
plot_cds.list <- cds.list

for(cell_line in names(plot_cds.list)){
  
  plot_cds.list[[cell_line]] <- plot_cds.list[[cell_line]][,!is.na(pData(plot_cds.list[[cell_line]])$Pseudotime_bin)]
  pData(plot_cds.list[[cell_line]])$Pseudotime_bin <- factor(pData(plot_cds.list[[cell_line]])$Pseudotime_bin,
                                                             levels = c("1","2","3","4","5","6","7","8","9","10"))
  
}


plot_cells(plot_cds.list[["A549"]], 
           reduction_method = "UMAP", color_cells_by = "Pseudotime_bin", label_cell_groups = FALSE,
           label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE, cell_size = 0.5, trajectory_graph_segment_size = 0.5) +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.25,"line"), 
        legend.key.height = unit(0.25,"line")) +
  guides(guides(colour = guide_legend(override.aes = list(size=1)))) +
  scale_color_manual("Pseudodose bin", values = c("1" = "#191970","2" = "#1C68CF","3" = "#238EC7","4" = "#2E8B57","5" = "#0FAE1D",
                                                  "6" = "#4FC300","7" = "#EEC900","8" = "#F99700","9" = "#FF5400","10" = "#FF0000")) +
  ggsave("supplemental_figure31_Aligned_A549_HDACi_UMAP_screen_by_Pseudotime.png", 
         dpi = 600, height = 1.5, width = 2.15)

plot_cells(plot_cds.list[["MCF7"]], reduction_method = "UMAP", color_cells_by = "Pseudotime_bin", label_cell_groups = FALSE,
           label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE, cell_size = 0.5, trajectory_graph_segment_size = 0.5) +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.25,"line"), 
        legend.key.height = unit(0.25,"line")) +
  guides(guides(colour = guide_legend(override.aes = list(size=1)))) +
  scale_color_manual("Pseudodose bin", values = c("1" = "#191970","2" = "#1C68CF","3" = "#238EC7","4" = "#2E8B57","5" = "#0FAE1D",
                                                  "6" = "#4FC300","7" = "#EEC900","8" = "#F99700","9" = "#FF5400","10" = "#FF0000")) +
  ggsave("supplemental_figure31_Aligned_MCF7_HDACi_UMAP_screen_by_Pseudotime.png", 
         dpi = 600, height = 1.5, width = 2.15)

# Create dataframes of cells by rescue and phenocopy experiment to determine the effect of modulation of Ac-CoA on HDACi Pseudotime
A549_rescue_df <- pData(cds.list[["A549"]]) %>%
  as.data.frame() %>% 
  filter(treatment %in% c("DMSO","Acetate","Citrate", "Pyruvate") & hdac != "Abexinostat") %>% 
  dplyr::select(hdac, hdac_dose, treatment, treatment_dose, Pseudotime, Pseudotime_bin) %>%
  filter(!is.na(Pseudotime_bin)) %>% 
  dplyr::mutate(Pseudotime_bin = as.character(Pseudotime_bin)) %>%
  group_by(Pseudotime_bin, hdac, hdac_dose, treatment, treatment_dose) %>% 
  dplyr::mutate(cells_per_bin = n(), combined_treatment = paste0(hdac,"_",hdac_dose,"_",treatment,"_",treatment_dose)) %>% 
  group_by(combined_treatment) %>%
  dplyr::mutate(fraction_of_cells_per_bin = cells_per_bin/n())

A549_rescue_df_overlapping_doses_only <- A549_rescue_df %>% filter(hdac_dose %in% c(0, 1000, 10000))

A549_rescue_df_subset <- A549_rescue_df_overlapping_doses_only %>% 
  dplyr::select(combined_treatment, Pseudotime_bin, fraction_of_cells_per_bin) %>%
  distinct()

A549_rescue_matrix <- dcast(A549_rescue_df_subset, combined_treatment~Pseudotime_bin)
row.names(A549_rescue_matrix) <- A549_rescue_matrix$combined_treatment
A549_rescue_matrix$combined_treatment <- NULL

A549_rescue_matrix[is.na(A549_rescue_matrix)] <- 0

sample_order <- c("DMSO_0_DMSO_0", "DMSO_0_Acetate_5", "DMSO_0_Acetate_10", "DMSO_0_Acetate_50", "DMSO_0_Acetate_100",
  "DMSO_0_Pyruvate_5", "DMSO_0_Pyruvate_10", "DMSO_0_Pyruvate_50", "DMSO_0_Pyruvate_100",
  "DMSO_0_Citrate_1", "DMSO_0_Citrate_2", "DMSO_0_Citrate_10","DMSO_0_Citrate_20",
  "Pracinostat_1000_DMSO_0", 
  "Pracinostat_1000_Acetate_5","Pracinostat_1000_Acetate_10","Pracinostat_1000_Acetate_50","Pracinostat_1000_Acetate_100",
  "Pracinostat_1000_Pyruvate_5","Pracinostat_1000_Pyruvate_10","Pracinostat_1000_Pyruvate_50","Pracinostat_1000_Pyruvate_100",
  "Pracinostat_1000_Citrate_1", "Pracinostat_1000_Citrate_2","Pracinostat_1000_Citrate_10","Pracinostat_1000_Citrate_20",
  "Pracinostat_10000_DMSO_0",
  "Pracinostat_10000_Acetate_5","Pracinostat_10000_Acetate_10","Pracinostat_10000_Acetate_50","Pracinostat_10000_Acetate_100",
  "Pracinostat_10000_Pyruvate_5","Pracinostat_10000_Pyruvate_10","Pracinostat_10000_Pyruvate_50","Pracinostat_10000_Pyruvate_100",
  "Pracinostat_10000_Citrate_1", "Pracinostat_10000_Citrate_2","Pracinostat_10000_Citrate_10","Pracinostat_10000_Citrate_20")

bin_order <- c("1","2","3","4","5","6","7","8","9","10")

A549_rescue_matrix <- A549_rescue_matrix[sample_order,bin_order]

pheatmap(A549_rescue_matrix, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "none",
         color = colorRampPalette(c("midnightblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1"))(25),
         file = "supplemental_figure31_A549_rescue_110519_key.png")

pheatmap(A549_rescue_matrix, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "none",
         color = colorRampPalette(c("midnightblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1"))(25),
         gaps_row = c(1,5,9,13,14,18,22,26,27,31,35,39),
         show_rownames = FALSE,
         show_colnames = FALSE,
         file = "supplemental_figure31_A549_rescue_110519.png",
         width = 2,
         height = 3)

A549_rescue_df %>% 
  ungroup() %>%
  filter(treatment_dose %in% c(0, 20, 100) & hdac_dose %in% c(0, 1000, 10000)) %>% 
  group_by(treatment) %>% 
  mutate(normalized_pseudotime = (Pseudotime + 1) - (mean(Pseudotime[hdac == "DMSO"]) + 1)) %>%
  filter(hdac != "DMSO") %>% 
  group_by(combined_treatment) %>%
  ungroup() %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("DMSO", "Acetate", 
                                                         "Pyruvate","Citrate"))) %>%
  ggplot() +
  geom_boxplot(aes(x = treatment, y = normalized_pseudotime), fill = "grey80", lwd = 0.1, outlier.shape = NA) +
  facet_wrap(~hdac_dose, ncol = 3, scales = "free_y") +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_x_discrete(
                    labels = c("DMSO" = "untreated","Acetate" = "100 mM acetate",
                               "Pyruvate" = "100 mM pyruvate","Citrate" =  "20 mM citrate")) +
  ylab("Normalized pseudodose") +
  monocle3::monocle_theme_opts() +
  ggsave("supplemental_figure31_A549_rescue_boxplots.png",dpi = 600, height = 1.75, width = 2.15)

##### Plot to examine q values
A549_rescue_df %>% 
  ungroup() %>%
  filter(treatment_dose %in% c(0, 20, 100) & hdac_dose %in% c(0, 1000, 10000)) %>% 
  group_by(treatment) %>% 
  mutate(normalized_pseudotime = (Pseudotime + 1) - (mean(Pseudotime[hdac == "DMSO"]) + 1)) %>%
  filter(hdac != "DMSO") %>% 
  ungroup() %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("DMSO", "Acetate", 
                                                         "Pyruvate","Citrate"))) %>%
  ggplot() +
  geom_boxplot(aes(x = treatment, y = normalized_pseudotime), fill = "grey80", lwd = 0.1, outlier.shape = NA) +
  facet_wrap(~hdac_dose, ncol = 3, scales = "free_y") +
  ggpubr::stat_compare_means(aes(x = treatment, y = normalized_pseudotime), method = "wilcox.test", ref.group = "DMSO", label = "p.signif") +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_x_discrete(#values = c("grey80","grey80","grey80","grey80"),
    labels = c("DMSO" = "untreated","Acetate" = "100 mM acetate",
               "Pyruvate" = "100 mM pyruvate","Citrate" =  "20 mM citrate")) +
  #ylim(NA, 15) +
  monocle3::monocle_theme_opts() +
  ggsave("supplemental_figure31_A549_rescue_boxplots_qvalue.png",height = 4, width = 10)


A549_phenocopy_df <- pData(cds.list[["A549"]]) %>%
  as.data.frame() %>% 
  filter(treatment %in% c("ACLY.inhibitor","ACSS2.inhibitor", "PDH.inhibitor") & hdac != "Abexinostat") %>% 
  dplyr::select(hdac, hdac_dose, treatment, treatment_dose, Pseudotime, Pseudotime_bin) %>%
  filter(!is.na(Pseudotime_bin)) %>% 
  dplyr::mutate(Pseudotime_bin = as.character(Pseudotime_bin), treatment = ifelse(treatment_dose == 0, "DMSO", treatment)) %>%
  group_by(Pseudotime_bin, hdac, hdac_dose, treatment, treatment_dose) %>% 
  dplyr::mutate(cells_per_bin = n(), combined_treatment = paste0(hdac,"_",hdac_dose,"_",treatment,"_",treatment_dose)) %>% 
  group_by(combined_treatment) %>%
  dplyr::mutate(fraction_of_cells_per_bin = cells_per_bin/n())

A549_phenocopy_df_overlapping_doses_only <- A549_phenocopy_df %>% filter(hdac_dose %in% c(0, 1000))

A549_phenocopy_df_subset <- A549_phenocopy_df_overlapping_doses_only %>% 
  dplyr::select(combined_treatment, Pseudotime_bin, fraction_of_cells_per_bin) %>%
  distinct()

A549_phenocopy_matrix <- dcast(A549_phenocopy_df_subset, combined_treatment~Pseudotime_bin)
row.names(A549_phenocopy_matrix) <- A549_phenocopy_matrix$combined_treatment
A549_phenocopy_matrix$combined_treatment <- NULL

A549_phenocopy_matrix[is.na(A549_phenocopy_matrix)] <- 0

sample_order <- c("DMSO_0_DMSO_0", 
                  "DMSO_0_ACLY.inhibitor_0.05", "DMSO_0_ACLY.inhibitor_0.25", "DMSO_0_ACLY.inhibitor_0.5", "DMSO_0_ACLY.inhibitor_2.5",
                  "DMSO_0_ACLY.inhibitor_5","DMSO_0_ACLY.inhibitor_25","DMSO_0_ACLY.inhibitor_50",
                  "DMSO_0_ACSS2.inhibitor_0.05", "DMSO_0_ACSS2.inhibitor_0.25", "DMSO_0_ACSS2.inhibitor_0.5", "DMSO_0_ACSS2.inhibitor_2.5",
                  "DMSO_0_ACSS2.inhibitor_5","DMSO_0_ACSS2.inhibitor_25","DMSO_0_ACSS2.inhibitor_50",
                  "DMSO_0_PDH.inhibitor_0.05", "DMSO_0_PDH.inhibitor_0.25", "DMSO_0_PDH.inhibitor_0.5", "DMSO_0_PDH.inhibitor_2.5",
                  "DMSO_0_PDH.inhibitor_5","DMSO_0_PDH.inhibitor_25","DMSO_0_PDH.inhibitor_50",
                  "Pracinostat_1000_DMSO_0", 
                  "Pracinostat_1000_ACLY.inhibitor_0.05", "Pracinostat_1000_ACLY.inhibitor_0.25", "Pracinostat_1000_ACLY.inhibitor_0.5", "Pracinostat_1000_ACLY.inhibitor_2.5",
                  "Pracinostat_1000_ACLY.inhibitor_5","Pracinostat_1000_ACLY.inhibitor_25","Pracinostat_1000_ACLY.inhibitor_50",
                  "Pracinostat_1000_ACSS2.inhibitor_0.05", "Pracinostat_1000_ACSS2.inhibitor_0.25", "Pracinostat_1000_ACSS2.inhibitor_0.5", "Pracinostat_1000_ACSS2.inhibitor_2.5",
                  "Pracinostat_1000_ACSS2.inhibitor_5","Pracinostat_1000_ACSS2.inhibitor_25","Pracinostat_1000_ACSS2.inhibitor_50",
                  "Pracinostat_1000_PDH.inhibitor_0.05", "Pracinostat_1000_PDH.inhibitor_0.25", "Pracinostat_1000_PDH.inhibitor_0.5", "Pracinostat_1000_PDH.inhibitor_2.5",
                  "Pracinostat_1000_PDH.inhibitor_5","Pracinostat_1000_PDH.inhibitor_25","Pracinostat_1000_PDH.inhibitor_50")

A549_phenocopy_matrix <- A549_phenocopy_matrix[sample_order,bin_order]

pheatmap(A549_phenocopy_matrix, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "none",
         color = colorRampPalette(c("midnightblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1"))(25),
         gaps_row = c(1,8,15,22,23,30,37,44),
         file = "supplemental_figure31_A549_Proportion_of_HDACi_and_inhibitor_exposed_cells_across_HDACi_trajectory_key.png")

pheatmap(A549_phenocopy_matrix, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "none",
         color = colorRampPalette(c("midnightblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1"))(25),
         show_rownames = FALSE,
         show_colnames = FALSE,
         gaps_row = c(1,8,15,22,23,30,37,44),
         file = "supplemental_figure31_A549_Proportion_of_HDACi_and_inhibitor_exposed_cells_across_HDACi_trajectory.png",
         width = 2,
         height = 3)


A549_phenocopy_df %>% 
  ungroup() %>%
  filter(treatment_dose %in% c(0,50)) %>% 
  group_by(treatment) %>% 
  mutate(normalized_pseudotime = (Pseudotime + 1) - (mean(Pseudotime[hdac == "DMSO"]) + 1)) %>%
  filter(hdac != "DMSO") %>% 
  group_by(combined_treatment) %>%
  ungroup() %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("DMSO", "ACLY.inhibitor", 
                                                         "ACSS2.inhibitor","PDH.inhibitor"))) %>%
  ggplot() +
  geom_boxplot(aes(x = treatment, y = normalized_pseudotime), fill = "grey80", lwd = 0.1, outlier.shape = NA) +
  facet_wrap(~hdac, scales = "free_y") +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_x_discrete(
    labels = c("DMSO" = "DMSO","ACLY.inhibitor" = "50 µM ACLYi",
               "ACSS2.inhibitor" = "50 µM ACSS2i","PDH.inhibitor" =  "50 µM PDHi")) +
  ylab("Normalized pseudodose") +
  monocle3::monocle_theme_opts() +
  ggsave("supplemental_figure31_A549_boxplots_pseudotime_distribution_pracinostat_inhibitors.png", 
         dpi = 600, height = 1.75, width = 1.2)

##### Plot to examine q-values
A549_phenocopy_df %>% 
  ungroup() %>%
  filter(treatment_dose %in% c(0,50)) %>% 
  group_by(treatment) %>% 
  mutate(normalized_pseudotime = (Pseudotime + 1) - (mean(Pseudotime[hdac == "DMSO"]) + 1)) %>%
  filter(hdac != "DMSO") %>% 
  ungroup() %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("DMSO", "ACLY.inhibitor", 
                                                         "ACSS2.inhibitor","PDH.inhibitor"))) %>%
  ggplot() +
  geom_boxplot(aes(x = treatment, y = normalized_pseudotime), fill = "grey80", lwd = 0.1, outlier.shape = NA) +
  facet_wrap(~hdac, scales = "free_y") +
  ggpubr::stat_compare_means(aes(x = treatment, y = normalized_pseudotime, label = "..p.adj.."), method = "wilcox.test", ref.group = "DMSO", label = "p.signif") +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_x_discrete(
    labels = c("DMSO" = "DMSO","ACLY.inhibitor" = "50 µM ACLYi",
               "ACSS2.inhibitor" = "50 µM ACSS2i","PDH.inhibitor" =  "50 µM PDHi")) +
  monocle3::monocle_theme_opts() +
  ggsave("supplemental_figure31_A549_boxplots_pseudotime_distribution_pracinostat_inhibitors_qvalue.png", 
         dpi = 600, height = 4, width = 7)

MCF7_rescue_df <- pData(cds.list[["MCF7"]]) %>%
  as.data.frame() %>% 
  filter(treatment %in% c("DMSO","Acetate","Citrate", "Pyruvate") & hdac != "Abexinostat") %>% 
  dplyr::select(hdac, hdac_dose, treatment, treatment_dose, Pseudotime, Pseudotime_bin) %>%
  filter(!is.na(Pseudotime_bin)) %>% 
  dplyr::mutate(Pseudotime_bin = as.character(Pseudotime_bin)) %>%
  group_by(Pseudotime_bin, hdac, hdac_dose, treatment, treatment_dose) %>% 
  dplyr::mutate(cells_per_bin = n(), combined_treatment = paste0(hdac,"_",hdac_dose,"_",treatment,"_",treatment_dose)) %>% 
  group_by(combined_treatment) %>%
  dplyr::mutate(fraction_of_cells_per_bin = cells_per_bin/n())

MCF7_rescue_df_overlapping_doses_only <- MCF7_rescue_df %>% filter(hdac_dose %in% c(0, 1000, 10000))

MCF7_rescue_df_subset <- MCF7_rescue_df_overlapping_doses_only %>% 
  dplyr::select(combined_treatment, Pseudotime_bin, fraction_of_cells_per_bin) %>%
  distinct()

MCF7_rescue_matrix <- dcast(MCF7_rescue_df_subset, combined_treatment~Pseudotime_bin)
row.names(MCF7_rescue_matrix) <- MCF7_rescue_matrix$combined_treatment
MCF7_rescue_matrix$combined_treatment <- NULL

MCF7_rescue_matrix[is.na(MCF7_rescue_matrix)] <- 0

sample_order <- c("DMSO_0_DMSO_0", "DMSO_0_Acetate_5", "DMSO_0_Acetate_10", "DMSO_0_Acetate_50", "DMSO_0_Acetate_100",
                  "DMSO_0_Pyruvate_5", "DMSO_0_Pyruvate_10", "DMSO_0_Pyruvate_50", "DMSO_0_Pyruvate_100",
                  "DMSO_0_Citrate_1", "DMSO_0_Citrate_2", "DMSO_0_Citrate_10","DMSO_0_Citrate_20",
                  "Pracinostat_1000_DMSO_0", 
                  "Pracinostat_1000_Acetate_5","Pracinostat_1000_Acetate_10","Pracinostat_1000_Acetate_50","Pracinostat_1000_Acetate_100",
                  "Pracinostat_1000_Pyruvate_5","Pracinostat_1000_Pyruvate_10","Pracinostat_1000_Pyruvate_50","Pracinostat_1000_Pyruvate_100",
                  "Pracinostat_1000_Citrate_1", "Pracinostat_1000_Citrate_2","Pracinostat_1000_Citrate_10","Pracinostat_1000_Citrate_20",
                  "Pracinostat_10000_DMSO_0",
                  "Pracinostat_10000_Acetate_5","Pracinostat_10000_Acetate_10","Pracinostat_10000_Acetate_50","Pracinostat_10000_Acetate_100",
                  "Pracinostat_10000_Pyruvate_5","Pracinostat_10000_Pyruvate_10","Pracinostat_10000_Pyruvate_50","Pracinostat_10000_Pyruvate_100",
                  "Pracinostat_10000_Citrate_1", "Pracinostat_10000_Citrate_2","Pracinostat_10000_Citrate_10","Pracinostat_10000_Citrate_20")

bin_order <- c("1","2","3","4","5","6","7","8","9","10")

MCF7_rescue_matrix <- MCF7_rescue_matrix[sample_order,bin_order]

pheatmap(MCF7_rescue_matrix, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "none",
         color = colorRampPalette(c("midnightblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1"))(25),
         file = "supplemental_figure31_MCF7_Proportion_of_HDACi_and_supplement_exposed_cells_across_HDACi_trajectory_Key.png")

pheatmap(MCF7_rescue_matrix, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "none",
         color = colorRampPalette(c("midnightblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1"))(25),
         gaps_row = c(1,5,9,13,14,18,22,26,27,31,35,39),
         show_rownames = FALSE,
         show_colnames = FALSE,
         file = "supplemental_figure31_MCF7_Proportion_of_HDACi_and_supplement_exposed_cells_across_HDACi_trajectory.png",
         width = 2,
         height = 3)

MCF7_rescue_df %>% 
  ungroup() %>%
  filter(treatment_dose %in% c(0, 20, 100) & hdac_dose %in% c(0, 1000, 10000)) %>% 
  group_by(treatment) %>% 
  mutate(normalized_pseudotime = (Pseudotime + 1) - (mean(Pseudotime[hdac == "DMSO"]) + 1)) %>%
  filter(hdac != "DMSO") %>% 
  group_by(combined_treatment) %>%
  ungroup() %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("DMSO", "Acetate", 
                                                         "Pyruvate","Citrate"))) %>%
  ggplot() +
  geom_boxplot(aes(x = treatment, y = normalized_pseudotime), fill = "grey80", lwd = 0.1, outlier.shape = NA) +
  facet_wrap(~hdac_dose, ncol = 3, scales = "free_y") +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_x_discrete(
    labels = c("DMSO" = "untreated","Acetate" = "100 mM acetate",
               "Pyruvate" = "100 mM pyruvate","Citrate" =  "20 mM citrate")) +
  ylab("Normalized pseudodose") +
  monocle3::monocle_theme_opts() +
  ggsave("supplemental_figure31_MCF7_boxplots_pseudotime_distribution_pracinostat_supplements.png", 
         dpi = 600, height = 1.75, width = 2.15)



#### Plot to examine q-values
MCF7_rescue_df %>% 
  ungroup() %>%
  filter(treatment_dose %in% c(0, 20, 100) & hdac_dose %in% c(0, 1000, 10000)) %>% 
  group_by(treatment) %>% 
  mutate(normalized_pseudotime = (Pseudotime + 1) - (mean(Pseudotime[hdac == "DMSO"]) + 1)) %>%
  filter(hdac != "DMSO") %>% 
  ungroup() %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("DMSO", "Acetate", 
                                                         "Pyruvate","Citrate"))) %>%
  ggplot() +
  geom_boxplot(aes(x = treatment, y = normalized_pseudotime), fill = "grey80", lwd = 0.1, outlier.shape = NA) +
  facet_wrap(~hdac_dose, ncol = 3, scales = "free_y") +
  ggpubr::stat_compare_means(aes(x = treatment, y = normalized_pseudotime, label = "..p.adj.."), method = "wilcox.test", ref.group = "DMSO", label = "p.signif") +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_x_discrete(
    labels = c("DMSO" = "untreated","Acetate" = "100 mM acetate",
               "Pyruvate" = "100 mM pyruvate","Citrate" =  "20 mM citrate")) +
  monocle3::monocle_theme_opts() +
  ggsave("supplemental_figure31_MCF7_boxplots_pseudotime_distribution_pracinostat_supplements.png", 
         dpi = 600, height = 4, width = 10)

MCF7_phenocopy_df <- pData(cds.list[["MCF7"]]) %>%
  as.data.frame() %>% 
  filter(treatment %in% c("ACLY.inhibitor","ACSS2.inhibitor", "PDH.inhibitor") & hdac != "Abexinostat") %>% 
  dplyr::select(hdac, hdac_dose, treatment, treatment_dose, Pseudotime, Pseudotime_bin) %>%
  filter(!is.na(Pseudotime_bin)) %>% 
  dplyr::mutate(Pseudotime_bin = as.character(Pseudotime_bin), treatment = ifelse(treatment_dose == 0, "DMSO", treatment)) %>%
  group_by(Pseudotime_bin, hdac, hdac_dose, treatment, treatment_dose) %>% 
  dplyr::mutate(cells_per_bin = n(), combined_treatment = paste0(hdac,"_",hdac_dose,"_",treatment,"_",treatment_dose)) %>% 
  group_by(combined_treatment) %>%
  dplyr::mutate(fraction_of_cells_per_bin = cells_per_bin/n())

MCF7_phenocopy_df_overlapping_doses_only <- MCF7_phenocopy_df %>% filter(hdac_dose %in% c(0, 1000))

MCF7_phenocopy_df_subset <- MCF7_phenocopy_df_overlapping_doses_only %>% 
  dplyr::select(combined_treatment, Pseudotime_bin, fraction_of_cells_per_bin) %>%
  distinct()

MCF7_phenocopy_matrix <- dcast(MCF7_phenocopy_df_subset, combined_treatment~Pseudotime_bin)
row.names(MCF7_phenocopy_matrix) <- MCF7_phenocopy_matrix$combined_treatment
MCF7_phenocopy_matrix$combined_treatment <- NULL

MCF7_phenocopy_matrix[is.na(MCF7_phenocopy_matrix)] <- 0

sample_order <- c("DMSO_0_DMSO_0", 
                  "DMSO_0_ACLY.inhibitor_0.05", "DMSO_0_ACLY.inhibitor_0.25", "DMSO_0_ACLY.inhibitor_0.5", "DMSO_0_ACLY.inhibitor_2.5",
                  "DMSO_0_ACLY.inhibitor_5","DMSO_0_ACLY.inhibitor_25","DMSO_0_ACLY.inhibitor_50",
                  "DMSO_0_ACSS2.inhibitor_0.05", "DMSO_0_ACSS2.inhibitor_0.25", "DMSO_0_ACSS2.inhibitor_0.5", "DMSO_0_ACSS2.inhibitor_2.5",
                  "DMSO_0_ACSS2.inhibitor_5","DMSO_0_ACSS2.inhibitor_25","DMSO_0_ACSS2.inhibitor_50",
                  "DMSO_0_PDH.inhibitor_0.05", "DMSO_0_PDH.inhibitor_0.25", "DMSO_0_PDH.inhibitor_0.5", "DMSO_0_PDH.inhibitor_2.5",
                  "DMSO_0_PDH.inhibitor_5","DMSO_0_PDH.inhibitor_25","DMSO_0_PDH.inhibitor_50",
                  "Pracinostat_1000_DMSO_0", 
                  "Pracinostat_1000_ACLY.inhibitor_0.05", "Pracinostat_1000_ACLY.inhibitor_0.25", "Pracinostat_1000_ACLY.inhibitor_0.5", "Pracinostat_1000_ACLY.inhibitor_2.5",
                  "Pracinostat_1000_ACLY.inhibitor_5","Pracinostat_1000_ACLY.inhibitor_25","Pracinostat_1000_ACLY.inhibitor_50",
                  "Pracinostat_1000_ACSS2.inhibitor_0.05", "Pracinostat_1000_ACSS2.inhibitor_0.25", "Pracinostat_1000_ACSS2.inhibitor_0.5", "Pracinostat_1000_ACSS2.inhibitor_2.5",
                  "Pracinostat_1000_ACSS2.inhibitor_5","Pracinostat_1000_ACSS2.inhibitor_25","Pracinostat_1000_ACSS2.inhibitor_50",
                  "Pracinostat_1000_PDH.inhibitor_0.05", "Pracinostat_1000_PDH.inhibitor_0.25", "Pracinostat_1000_PDH.inhibitor_0.5", "Pracinostat_1000_PDH.inhibitor_2.5",
                  "Pracinostat_1000_PDH.inhibitor_5","Pracinostat_1000_PDH.inhibitor_25","Pracinostat_1000_PDH.inhibitor_50")

MCF7_phenocopy_matrix <- MCF7_phenocopy_matrix[sample_order,bin_order]

pheatmap(MCF7_phenocopy_matrix, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "none",
         color = colorRampPalette(c("midnightblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1"))(25),
         gaps_row = c(1,8,15,22,23,30,37,44),
         file = "supplemental_figure31_MCF7_Proportion_of_HDACi_and_inhibitor_exposed_cells_across_HDACi_trajectory_key.png")

pheatmap(MCF7_phenocopy_matrix, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "none",
         color = colorRampPalette(c("midnightblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1"))(25),
         show_rownames = FALSE,
         show_colnames = FALSE,
         gaps_row = c(1,8,15,22,23,30,37,44),
         file = "supplemental_figure31_MCF7_Proportion_of_HDACi_and_inhibitor_exposed_cells_across_HDACi_trajectory.png",
         width = 2,
         height = 3)


MCF7_phenocopy_df %>% 
  ungroup() %>%
  filter(treatment_dose %in% c(0,50)) %>% 
  group_by(treatment) %>% 
  mutate(normalized_pseudotime = (Pseudotime + 1) - (mean(Pseudotime[hdac == "DMSO"]) + 1)) %>%
  filter(hdac != "DMSO") %>% 
  group_by(combined_treatment) %>%
  ungroup() %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("DMSO", "ACLY.inhibitor", 
                                                         "ACSS2.inhibitor","PDH.inhibitor"))) %>%
  ggplot() +
  geom_boxplot(aes(x = treatment, y = normalized_pseudotime), fill = "grey80", lwd = 0.1, outlier.shape = NA) +
  facet_wrap(~hdac, scales = "free_y") +

  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_x_discrete(#values = c("grey80","grey80","grey80","grey80"),
    labels = c("DMSO" = "DMSO","ACLY.inhibitor" = "50 µM ACLYi",
               "ACSS2.inhibitor" = "50 µM ACSS2i","PDH.inhibitor" =  "50 µM PDHi")) +
  #ylim(0, 15) +
  ylab("Normalized pseudodose") +
  monocle3::monocle_theme_opts() +
  ggsave("supplemental_figure31_MCF7_boxplots_pseudotime_distribution_pracinostat_inhibitors.png", 
         dpi = 600, height = 1.75, width = 1.2)

MCF7_phenocopy_df %>% 
  ungroup() %>%
  filter(treatment_dose %in% c(0,50)) %>% 
  group_by(treatment) %>% 
  mutate(normalized_pseudotime = (Pseudotime + 1) - (mean(Pseudotime[hdac == "DMSO"]) + 1)) %>%
  filter(hdac != "DMSO") %>% 
  ungroup() %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("DMSO", "ACLY.inhibitor", 
                                                         "ACSS2.inhibitor","PDH.inhibitor"))) %>%
  ggplot() +
  geom_boxplot(aes(x = treatment, y = normalized_pseudotime), fill = "grey80", lwd = 0.1, outlier.shape = NA) +
  facet_wrap(~hdac, scales = "free_y") +
  ggpubr::stat_compare_means(aes(x = treatment, y = normalized_pseudotime, label = "..p.adj.."), method = "wilcox.test", ref.group = "DMSO", label = "p.signif") +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_x_discrete(#values = c("grey80","grey80","grey80","grey80"),
    labels = c("DMSO" = "DMSO","ACLY.inhibitor" = "50 µM ACLYi",
               "ACSS2.inhibitor" = "50 µM ACSS2i","PDH.inhibitor" =  "50 µM PDHi")) +
  #ylim(0, 15) +
  monocle3::monocle_theme_opts() +
  ggsave("supplemental_figure31_MCF7_boxplots_pseudotime_distribution_pracinostat_inhibitors.png", 
         dpi = 600, height = 4, width = 10)

