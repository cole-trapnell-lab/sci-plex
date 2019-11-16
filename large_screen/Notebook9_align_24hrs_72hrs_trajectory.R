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
  library(plyr)
  library(dplyr)
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

# Set DelayedArray Parameters
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size = 1000e7)

# User defined functions

setwd(paste0(path_to_github,
             "large_screen",
             sep = ""))

# Path to cds containing all the data
path_to_cds = "cds.RDS"

# Read in CDS and extract cells passing filters and human transcripts
cds = readRDS(path_to_cds)

append_umap_coordinates = function(cds) {
  num_dim = dim(cds@reducedDims[["UMAP"]])[2]
  for (i in seq(1, num_dim)) {
    new_col_name = paste("umap", i, sep = "")
    colData(cds)[[new_col_name]] = cds@reducedDims[["UMAP"]][, i]
  }
  return(cds)
}


# Genes in the human transcriptome
human.genes = grepl("^ENSG", rowData(cds)$id)


# Cells that passed hash filters
validly_labeled_cells =
  colData(cds) %>%
  as.data.frame() %>%
  filter(!is.na(top_oligo_P),!is.na(top_oligo_W)) %>%
  filter(
    hash_umis_P >= 5,
    hash_umis_W >= 5,
    top_to_second_best_ratio_P >= 5,
    top_to_second_best_ratio_W >= 5
  ) %>%
  filter(rt_well != 1) %>%
  pull(cell)

cds = cds[human.genes, validly_labeled_cells]


compounds_profiled_at_72hrs =
  colData(cds) %>%
  as.data.frame() %>%
  filter(time_point == 72) %>%
  group_by(product_name) %>%
  pull(product_name) %>%
  unique()

A549_cells =
  colData(cds) %>%
  as.data.frame() %>%
  filter(!is.na(cell_type),
         cell_type == "A549") %>%
  pull(cell)

sciPlex_A549_cds <- cds[, A549_cells]
sciPlex_A549_cds <-
  sciPlex_A549_cds[, colData(sciPlex_A549_cds)$product_name %in% compounds_profiled_at_72hrs]

sciPlex_A549_24_72_hr_cds <- sciPlex_A549_cds

trajectory_drugs <-
  c(
    "S0000",
    "S2693",
    "S1045",
    "S1096",
    "S1053",
    "S1030",
    "S2170",
    "S2759",
    "S1085",
    "S1095",
    "S1194",
    "S2818",
    "S2244",
    "S1090",
    "S1515",
    "S1122",
    "S2779",
    "S8567"
  )

trajectory_drugs <-
  unique(pData(sciPlex_A549_24_72_hr_cds)[pData(sciPlex_A549_24_72_hr_cds)$catalog_number %in% trajectory_drugs, ]$product_name)

drugs_72_hr <- (
  pData(sciPlex_A549_24_72_hr_cds) %>%
    as.data.frame() %>%
    filter(time_point == 72) %>%
    group_by(product_name) %>%
    summarize(n = n())
)$product_name

HDAC_trajectory_drugs <-
  drugs_72_hr[drugs_72_hr %in% trajectory_drugs]

sciPlex_24_72hr_HDAC_subset <- sciPlex_A549_24_72_hr_cds
sciPlex_24_72hr_HDAC_subset <-
  sciPlex_24_72hr_HDAC_subset[, pData(sciPlex_24_72hr_HDAC_subset)$product_name %in% HDAC_trajectory_drugs]

pData(sciPlex_24_72hr_HDAC_subset)$viability <- NULL

cells_by_time_point <-
  pData(sciPlex_24_72hr_HDAC_subset) %>%
  as.data.frame() %>%
  group_by(time_point, product_name) %>%
  summarize(n = n())

Abexinostat_24hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 24 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Abexinostat (PCI-24781)", ])
Belinostat_24hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 24 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Belinostat (PXD101)", ])
Dacinostat_24hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 24 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Dacinostat (LAQ824)", ])
Mocetinostat_24hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 24 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Mocetinostat (MGCD0103)", ])
Panobinostat_24hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 24 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Panobinostat (LBH589)", ])
Pracinostat_24hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 24 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Pracinostat (SB939)", ])
Quisinostat_24hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 24 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Quisinostat (JNJ-26481585) 2HCl", ])
Tucidinostat_24hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 24 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Tucidinostat (Chidamide)", ])
Vehicle_24hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 24 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Vehicle", ])

Abexinostat_72hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 72 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Abexinostat (PCI-24781)", ])
Belinostat_72hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 72 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Belinostat (PXD101)", ])
Dacinostat_72hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 72 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Dacinostat (LAQ824)", ])
Mocetinostat_72hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 72 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Mocetinostat (MGCD0103)", ])
Panobinostat_72hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 72 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Panobinostat (LBH589)", ])
Pracinostat_72hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 72 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Pracinostat (SB939)", ])
Quisinostat_72hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 72 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Quisinostat (JNJ-26481585) 2HCl", ])
Tucidinostat_72hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 72 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Tucidinostat (Chidamide)", ])
Vehicle_72hr_cells <-
  row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)$time_point == 72 &
                                                 pData(sciPlex_24_72hr_HDAC_subset)$product_name == "Vehicle", ])


set.seed(1442)
downsampled_Abexinostat_cells_at_72_hrs <-
  sample(Abexinostat_72hr_cells, cells_by_time_point[cells_by_time_point$time_point == 24
                                                     &
                                                       cells_by_time_point$product_name == "Abexinostat (PCI-24781)", ]$n)
Abexinostat_cells <-
  c(Abexinostat_24hr_cells,
    downsampled_Abexinostat_cells_at_72_hrs)

set.seed(1442)
downsampled_Belinostat_cells_at_72_hrs <-
  sample(Belinostat_72hr_cells, cells_by_time_point[cells_by_time_point$time_point == 24
                                                    &
                                                      cells_by_time_point$product_name == "Belinostat (PXD101)", ]$n)
Belinostat_cells <-
  c(Belinostat_24hr_cells,
    downsampled_Belinostat_cells_at_72_hrs)

set.seed(1442)
downsampled_Dacinostat_cells_at_72_hrs <-
  sample(Dacinostat_72hr_cells, cells_by_time_point[cells_by_time_point$time_point == 24
                                                    &
                                                      cells_by_time_point$product_name == "Dacinostat (LAQ824)", ]$n)
Dacinostat_cells <-
  c(Dacinostat_24hr_cells,
    downsampled_Dacinostat_cells_at_72_hrs)

set.seed(1442)
downsampled_Mocetinostat_cells_at_72_hrs <-
  sample(Mocetinostat_72hr_cells, cells_by_time_point[cells_by_time_point$time_point == 24
                                                      &
                                                        cells_by_time_point$product_name == "Mocetinostat (MGCD0103)", ]$n)
Mocetinostat_cells <-
  c(Mocetinostat_24hr_cells,
    downsampled_Mocetinostat_cells_at_72_hrs)

set.seed(1442)
downsampled_Panobinostat_cells_at_24_hrs <-
  sample(Panobinostat_24hr_cells, cells_by_time_point[cells_by_time_point$time_point == 72
                                                      &
                                                        cells_by_time_point$product_name == "Panobinostat (LBH589)", ]$n)
Panobinostat_cells <-
  c(Panobinostat_72hr_cells,
    downsampled_Panobinostat_cells_at_24_hrs)

set.seed(1442)
downsampled_Pracinostat_cells_at_72_hrs <-
  sample(Pracinostat_72hr_cells, cells_by_time_point[cells_by_time_point$time_point == 24
                                                     &
                                                       cells_by_time_point$product_name == "Pracinostat (SB939)", ]$n)
Pracinostat_cells <-
  c(Pracinostat_24hr_cells,
    downsampled_Pracinostat_cells_at_72_hrs)

set.seed(1442)
downsampled_Quisinostat_cells_at_24_hrs <-
  sample(Quisinostat_24hr_cells, cells_by_time_point[cells_by_time_point$time_point == 72
                                                     &
                                                       cells_by_time_point$product_name == "Quisinostat (JNJ-26481585) 2HCl", ]$n)
Quisinostat_cells <-
  c(Quisinostat_72hr_cells,
    downsampled_Quisinostat_cells_at_24_hrs)

set.seed(1442)
downsampled_Tucidinostat_cells_at_72_hrs <-
  sample(Tucidinostat_72hr_cells, cells_by_time_point[cells_by_time_point$time_point == 24
                                                      &
                                                        cells_by_time_point$product_name == "Tucidinostat (Chidamide)", ]$n)
Tucidinostat_cells <-
  c(Tucidinostat_24hr_cells,
    downsampled_Tucidinostat_cells_at_72_hrs)

set.seed(1442)
downsampled_Vehicle_cells_at_24_hrs <-
  sample(Vehicle_24hr_cells, cells_by_time_point[cells_by_time_point$time_point == 72
                                                 &
                                                   cells_by_time_point$product_name == "Vehicle", ]$n)
Vehicle_cells <-
  c(Vehicle_72hr_cells, downsampled_Vehicle_cells_at_24_hrs)


sciPlex_24_72hr_HDAC_subset <-
  sciPlex_24_72hr_HDAC_subset[, c(
    Abexinostat_cells,
    Belinostat_cells,
    Mocetinostat_cells,
    Dacinostat_cells,
    Panobinostat_cells,
    Pracinostat_cells,
    Quisinostat_cells,
    Tucidinostat_cells,
    Vehicle_cells
  )]

sciPlex_24_72hr_HDAC_subset <-
  estimate_size_factors(sciPlex_24_72hr_HDAC_subset)

sciPlex_24_HDAC_subset_counts = count_cells(sciPlex_24_72hr_HDAC_subset[, pData(sciPlex_24_72hr_HDAC_subset)$time_point == 24])
sciPlex_24_HDAC_subset_models = fit_dose_response_models(sciPlex_24_HDAC_subset_counts)
sciPlex_24_HDAC_subset_model_stats = extract_dose_model_stats(sciPlex_24_HDAC_subset_counts,
                                                              sciPlex_24_HDAC_subset_models)

sciPlex_72_HDAC_subset_counts = count_cells(sciPlex_24_72hr_HDAC_subset[, pData(sciPlex_24_72hr_HDAC_subset)$time_point == 72])
sciPlex_72_HDAC_subset_models = fit_dose_response_models(sciPlex_72_HDAC_subset_counts)
sciPlex_72_HDAC_subset_model_stats = extract_dose_model_stats(sciPlex_72_HDAC_subset_counts,
                                                              sciPlex_72_HDAC_subset_models)

estimate_viability_pData <- function (cds, sci_plex_models)
{
  viability_df =
    pData(cds) %>%
    as.data.frame() %>%
    group_by(cell_type, sample, treatment) %>%
    do(estimate_viability_helper(., sci_plex_models))
  
  row.names(viability_df) = viability_df$cell
  pData(cds)$viability = pmin(1, viability_df[pData(cds)$cell,]$viability)
  return(pData(cds))
}

sciPlex_A549_24hr_viability_pData <-
  estimate_viability_pData(sciPlex_24_72hr_HDAC_subset[, pData(sciPlex_24_72hr_HDAC_subset)$time_point == 24],
                           sciPlex_24_HDAC_subset_models)

sciPlex_A549_72hr_viability_pData <-
  estimate_viability_pData(sciPlex_24_72hr_HDAC_subset[, pData(sciPlex_24_72hr_HDAC_subset)$time_point == 72],
                           sciPlex_72_HDAC_subset_models)

combined_viability_pData <- rbind(sciPlex_A549_24hr_viability_pData,
                                  sciPlex_A549_72hr_viability_pData)

viability_column =
  left_join(
    colData(sciPlex_24_72hr_HDAC_subset) %>%
      as.data.frame(),
    combined_viability_pData[, c("cell", "viability")] %>%
      as.data.frame(),
    by = "cell"
  ) %>%
  dplyr::select(cell, viability)

if (identical(colnames(sciPlex_24_72hr_HDAC_subset) ,
              viability_column$cell)) {
  colData(sciPlex_24_72hr_HDAC_subset)$viability = viability_column$viability
} else{
  message("Make sure that rownames match before joining")
}


pData(sciPlex_24_72hr_HDAC_subset)$viability[is.na(pData(sciPlex_24_72hr_HDAC_subset)$viability)] <-
  1

sciPlex_24_72hr_HDAC_subset <-
  estimate_cell_cycle(
    sciPlex_24_72hr_HDAC_subset,
    g1s_markers = cc.genes$s.genes,
    g2m_markers = cc.genes$g2m.genes
  )


# After regresing out the effects of number of UMIs and replicate, perform PCA and return 25 PCs
sciPlex_24_72hr_HDAC_subset <-
  preprocess_cds(
    sciPlex_24_72hr_HDAC_subset,
    method = 'PCA',
    num_dim = 25,
    norm_method = 'log',
    residual_model_formula_str = "~log(n.umi) + replicate + viability + proliferation_index",
    verbose = T,
    cores = 20
  )



cells_list <- list()
time_point <-
  unique(colData(sciPlex_24_72hr_HDAC_subset)[, "time_point"])

for (t in time_point) {
  cells_list[[as.character(t)]] <-
    row.names(pData(sciPlex_24_72hr_HDAC_subset)[pData(sciPlex_24_72hr_HDAC_subset)[, "time_point"] == t, ])
}

PCA_cells_1 = sciPlex_24_72hr_HDAC_subset@reducedDims[["PCA"]][cells_list[[1]], ]
PCA_cells_2 = sciPlex_24_72hr_HDAC_subset@reducedDims[["PCA"]][cells_list[[2]], ]



corrected_PCA <- scran::mnnCorrect(t(PCA_cells_1),
                                   t(PCA_cells_2))

corrected_PCA <- t(do.call("cbind", corrected_PCA[["corrected"]]))

row.names(corrected_PCA) <- c(row.names(PCA_cells_1),
                              row.names(PCA_cells_2))

corrected_PCA <-
  corrected_PCA[colnames(sciPlex_24_72hr_HDAC_subset), ]
sciPlex_24_72hr_HDAC_subset@reducedDims[["PCA"]] = as.matrix(corrected_PCA)




# Perform UMAP dimensionality reduction
sciPlex_24_72hr_HDAC_subset <-
  reduce_dimension(
    sciPlex_24_72hr_HDAC_subset,
    max_components = 2,
    reduction_method = 'UMAP',
    umap.metric = 'cosine',
    umap.n_neighbors = 10,
    umap.min_dist = 0.1,
    verbose = TRUE,
    cores = 20
  )


sciPlex_24_72hr_HDAC_subset <-
  cluster_cells(sciPlex_24_72hr_HDAC_subset, reduction_method = "PCA")

colData(sciPlex_24_72hr_HDAC_subset)$Cluster = clusters(sciPlex_24_72hr_HDAC_subset, reduction_method =
                                                          "PCA")

sciPlex_24_72hr_HDAC_subset <-
  cluster_cells(sciPlex_24_72hr_HDAC_subset, reduction_method = "UMAP")
colData(sciPlex_24_72hr_HDAC_subset)$louvain_component = sciPlex_24_72hr_HDAC_subset@clusters[["UMAP"]]$partitions

graph_parameters = list()
graph_parameters[["minimal_branch_len"]] = 10
graph_parameters[["ncenter"]] = 250

sciPlex_24_72hr_HDAC_subset <-
  learn_graph(
    sciPlex_24_72hr_HDAC_subset,
    learn_graph_control = graph_parameters,
    use_partition = T,
    close_loop = F
  )
sciPlex_24_72hr_HDAC_subset = append_umap_coordinates(sciPlex_24_72hr_HDAC_subset)



ggplot(
  pData(sciPlex_24_72hr_HDAC_subset) %>% as.data.frame(),
  aes(
    x = umap1,
    y = umap2,
    color = as.factor(time_point)
  )
) +
  geom_point(size = 0.1, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line")
  ) +
  scale_color_manual("Time point", values = c("24" = "brown4", "72" = "navy")) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure26_sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_alignment_HDAC_subset.png",
    width = 2,
    height = 1.5,
    dpi = 600
  )






################ new

# Get the closest vertex for every cell
pData(sciPlex_24_72hr_HDAC_subset)$closest_vertex <-
  sciPlex_24_72hr_HDAC_subset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[, 1]

ordering_summary =
  pData(sciPlex_24_72hr_HDAC_subset) %>%
  as.data.frame() %>%
  dplyr::group_by(closest_vertex) %>%
  dplyr::count(dose) %>%
  dplyr::mutate(total_cells = sum(n),
                fraction_dose = n / total_cells)

root_nodes =
  ordering_summary %>%
  filter(dose == 0 & fraction_dose > 0.5)

root_nodes = root_nodes$closest_vertex

pData(sciPlex_24_72hr_HDAC_subset)$root_node = pData(sciPlex_24_72hr_HDAC_subset)$closest_vertex %in% root_nodes

root_cells  =
  colData(sciPlex_24_72hr_HDAC_subset) %>%
  as.data.frame() %>%
  filter(root_node) %>%
  pull(cell) %>%
  as.character()

sciPlex_24_72hr_HDAC_subset <-
  order_cells(sciPlex_24_72hr_HDAC_subset, root_cells = root_cells)

pData(sciPlex_24_72hr_HDAC_subset)$Pseudotime = sciPlex_24_72hr_HDAC_subset@principal_graph_aux[["UMAP"]]$pseudotime

pData(sciPlex_24_72hr_HDAC_subset)$Pseudotime_bin =
  cut(
    pData(sciPlex_24_72hr_HDAC_subset)$Pseudotime,
    breaks = quantile(pData(sciPlex_24_72hr_HDAC_subset)$Pseudotime, seq(0, 1, 0.1)),
    labels = F
  )



vehicle.df =
  pData(sciPlex_24_72hr_HDAC_subset) %>%
  as.data.frame() %>%
  dplyr::select(Pseudotime,
                dose_character,
                product_name,
                cell_type,
                time_point) %>%
  filter(product_name == "Vehicle") %>%
  dplyr::select(-product_name)

non.vehicle.df =
  pData(sciPlex_24_72hr_HDAC_subset)  %>%
  as.data.frame() %>%
  dplyr::select(Pseudotime,
                dose_character,
                product_name,
                cell_type,
                time_point) %>%
  filter(product_name != "Vehicle")

unique.keys.df =
  pData(sciPlex_24_72hr_HDAC_subset)  %>%
  as.data.frame() %>%
  dplyr::select(cell_type, product_name) %>%
  unique()

duplicated.vehicle.df = inner_join(unique.keys.df, vehicle.df, by = "cell_type")

duplicated.vehicle.df =
  duplicated.vehicle.df %>%
  filter(product_name != "Vehicle") %>%
  dplyr::select(cell_type,
                product_name,
                Pseudotime,
                dose_character,
                time_point)

non.vehicle.df =
  non.vehicle.df %>%
  dplyr::select(cell_type,
                product_name,
                Pseudotime,
                dose_character,
                time_point)

combined.df = rbind(duplicated.vehicle.df, non.vehicle.df)



combined.df %>%
  filter(is.finite(Pseudotime)) %>%
  ggplot() +
  ggridges::geom_density_ridges(
    data = combined.df[combined.df$time_point == 24, ],
    aes(x = Pseudotime, y = dose_character),
    linetype = "dashed",
    # fill =  "grey80",
    size = 0.1,
    alpha = 0.4,
    scale = 1.5,
    color = "gray80"
  ) + # used to be 0.5
  ggridges::geom_density_ridges(
    data = combined.df[combined.df$time_point == 72, ],
    aes(x = Pseudotime, y = dose_character, fill = dose_character),
    # fill = "black",
    size = 0.1,
    alpha = 0.8,
    scale = 1.5
  ) + # used to be 1
  monocle:::monocle_theme_opts() +
  facet_wrap( ~ product_name, ncol = 4) +
  scale_fill_manual(
    values = c(
      "0" = "gray",
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0.5, "lines")
  ) +
  ggsave(
    "supplemental_figure26__24-72hr_ridge_plot_2.png",
    height = 1.5 ,
    width = 2.25,
    unit = "in",
    dpi = 1200
  )



combined.df %>%
  filter(is.finite(Pseudotime)) %>%
  ggplot() +
  ggridges::geom_density_ridges(aes(x = Pseudotime, y = dose_character, fill = dose_character),
                                size = 0.2) +
  monocle:::monocle_theme_opts() +
  facet_wrap(product_name ~ time_point, ncol = 8) +
  scale_fill_manual(
    values = c(
      "0" = "gray",
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ggsave(
    "Figure4_24-72hr_ridge_plot_w_titles.pdf",
    height = 3 ,
    width = 8,
    unit = "in"
  )


pData(sciPlex_24_72hr_HDAC_subset)$Pseudotime_bin <-
  as.character(pData(sciPlex_24_72hr_HDAC_subset)$Pseudotime_bin)
pData(sciPlex_24_72hr_HDAC_subset)$Pseudotime_bin <-
  factor(
    pData(sciPlex_24_72hr_HDAC_subset)$Pseudotime_bin,
    levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
  )



ggplot(
  pData(sciPlex_24_72hr_HDAC_subset) %>% as.data.frame(),
  aes(x = umap1, y = umap2, color = as.factor(dose))
) +
  geom_point(size = 0.2, stroke = 0) +
  facet_wrap( ~ time_point) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    axis.text = element_text(size = 4),
    legend.key.width = unit(0.5, "line"),
    legend.key.height = unit(0.4, "line"),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_color_manual(
    "Dose [nM]",
    values = c(
      "0" = "gray",
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1.5)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure26_A549_24_72_hr_aligned_HDAC_trajectory_by_Pseudotime.png",
    width = 2,
    height = 1.5,
    dpi = 600
  )



ica_space_df <-
  t(sciPlex_24_72hr_HDAC_subset@principal_graph_aux[["UMAP"]]$dp_mst) %>%
  as.data.frame() %>% dplyr::select(prin_graph_dim_1 = V1,
                                    prin_graph_dim_2 = V2) %>%
  dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))

edge_df <-
  sciPlex_24_72hr_HDAC_subset@principal_graph[["UMAP"]] %>%
  igraph::as_data_frame() %>%
  dplyr::select_(source = "from", target = "to") %>%
  dplyr::left_join(
    ica_space_df %>%
      dplyr::select(
        source = sample_name,
        source_prin_graph_dim_1 = prin_graph_dim_1,
        source_prin_graph_dim_2 = prin_graph_dim_2
      ),
    by = "source"
  ) %>%
  dplyr::left_join(
    ica_space_df %>%
      dplyr::select(
        target = sample_name,
        target_prin_graph_dim_1 = prin_graph_dim_1,
        target_prin_graph_dim_2 = prin_graph_dim_2
      ),
    by = "target"
  )


colData(sciPlex_24_72hr_HDAC_subset)$product_name_alt <-
  plyr::revalue((pData(sciPlex_24_72hr_HDAC_subset)$product_name),
                c(
                  "Abexinostat (PCI-24781)" = 'Abexinostat',
                  "Pracinostat (SB939)" = 'Pracinostat',
                  "Mocetinostat (MGCD0103)" = 'Mocetinostat',
                  "Tucidinostat (Chidamide)" = 'Tucidinostat',
                  "AR-42" = 'AR-42',
                  "Tacedinaline (CI994)" = 'Tacedinaline',
                  "CUDC-101" = 'CUDC-101',
                  "M344" = 'M344',
                  "CUDC-907" = 'CUDC-907',
                  "Dacinostat (LAQ824)" = 'Dacinostat',
                  "Belinostat (PXD101)" = 'Belinostat',
                  "vehicle" = 'Vehicle',
                  "Quisinostat (JNJ-26481585) 2HCl" = 'Quisinostat',
                  "Trichostatin A (TSA)" = 'TSA',
                  "Resminostat" = 'Resminostat',
                  "Entinostat (MS-275)" = 'Entinostat',
                  "Givinostat (ITF2357)" = 'Givinostat',
                  "Panobinostat (LBH589)" = 'Panobinostat',
                  "Vehicle" = 'Vehicle'
                )
  )

ggplot() +
  geom_point(
    data =
      colData(sciPlex_24_72hr_HDAC_subset) %>%
      as.data.frame() %>%
      filter(vehicle) %>%
      dplyr::select(-product_name_alt),
    aes(x = umap1, y = umap2),
    color = "grey80",
    size = .5,
    stroke = 0
  ) +
  geom_point(
    data =
      colData(sciPlex_24_72hr_HDAC_subset) %>%
      as.data.frame() %>%
      filter(!vehicle),
    aes(
      x = umap1,
      y = umap2,
      color = as.factor(time_point)
    ),
    size = .5,
    stroke = 0
  ) +
  geom_segment(
    data = edge_df,
    aes_string(
      x = "source_prin_graph_dim_1",
      y = "source_prin_graph_dim_2",
      xend = "target_prin_graph_dim_1",
      yend = "target_prin_graph_dim_2"
    ),
    size = .25,
    linetype = "solid",
    na.rm = TRUE
  ) +
  facet_wrap( ~ product_name_alt, nrow = 2) +
  scale_color_brewer("Time\nPoint", palette = "Set1") +
  monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  ggsave(
    "supplemental_figure26_UMAP_time_point_trajectory.png",
    width = 5,
    height = 3,
    unit = "in",
    dpi = 900
  )


ggplot() +
  geom_point(
    data =
      colData(sciPlex_24_72hr_HDAC_subset) %>%
      as.data.frame() %>%
      filter(is.finite(Pseudotime)),
    aes(
      x = umap1,
      y = umap2,
      color = as.numeric(Pseudotime)
    ),
    size = .5,
    stroke = 0
  ) +
  geom_segment(
    data = edge_df,
    aes_string(
      x = "source_prin_graph_dim_1",
      y = "source_prin_graph_dim_2",
      xend = "target_prin_graph_dim_1",
      yend = "target_prin_graph_dim_2"
    ),
    size = .25,
    linetype = "solid",
    na.rm = TRUE
  ) +
  facet_wrap( ~ time_point, nrow = 1) +
  scale_color_viridis_c("Pseudodose") +
  monocle_theme_opts() +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    strip.text.x = element_blank()
  ) +
  guides(color = guide_colourbar(barwidth = 0.75, barheight = 2)) +
  xlab("Component 1") +
  ylab("Component 2") +
  ggsave(
    "supplemental_figure26_UMAP_time_point_trajectory_psedodose.png",
    width = 4,
    height = 2,
    unit = "in",
    dpi = 900
  )
