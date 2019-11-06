# devtools::install_github(repo = 'cole-trapnell-lab/monocle3', ref = "d4a9a35")

# Input path to github directory
path_to_github = "~/sci-Plex/"
path_to_monocle3 = paste0(path_to_github,
                          "monocle3_d4a9a35/monocle3/",
                          sep = "")

# Set directory for sciPlex bin
bin_directory = paste0(path_to_github,
                       "bin",
                       sep = "")
suppressPackageStartupMessages({
  library(devtools)
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
  load_all(path_to_monocle3)
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

dim(sciPlex_A549_cds)

sciPlex_A549_24hr_counts = count_cells(sciPlex_A549_cds[, pData(sciPlex_A549_cds)$time_point == 24])
sciPlex_A549_24hr_models = fit_dose_response_models(sciPlex_A549_24hr_counts)
sciPlex_A549_24hr_model_stats = extract_dose_model_stats(sciPlex_A549_24hr_counts, sciPlex_A549_24hr_models)

sciPlex_A549_72hr_counts = count_cells(sciPlex_A549_cds[, pData(sciPlex_A549_cds)$time_point == 72])
sciPlex_A549_72hr_models = fit_dose_response_models(sciPlex_A549_72hr_counts)
sciPlex_A549_72hr_model_stats = extract_dose_model_stats(sciPlex_A549_72hr_counts, sciPlex_A549_72hr_models)



estimate_viability_pData <- function (cds, sci_Plex_models)
{
  viability_df =
    colData(cds) %>%
    as.data.frame() %>%
    group_by(cell_type, sample, treatment) %>%
    do(estimate_viability_helper(., sci_Plex_models))
  row.names(viability_df) = viability_df$cell
  pData(cds)$viability = pmin(1, viability_df[pData(cds)$cell, ]$viability)
  return(pData(cds))
}

sciPlex_A549_24hr_viability_pData <-
  estimate_viability_pData(sciPlex_A549_cds[, pData(sciPlex_A549_cds)$time_point == 24],
                           sciPlex_A549_24hr_models)

sciPlex_A549_72hr_viability_pData <-
  estimate_viability_pData(sciPlex_A549_cds[, pData(sciPlex_A549_cds)$time_point == 72],
                           sciPlex_A549_72hr_models)


combined_viability_pData <- rbind(sciPlex_A549_24hr_viability_pData,
                                  sciPlex_A549_72hr_viability_pData)

viability_column =
  left_join(
    colData(sciPlex_A549_cds) %>%
      as.data.frame(),
    combined_viability_pData[, c("cell", "viability")] %>%
      as.data.frame(),
    by = "cell"
  ) %>%
  dplyr::select(cell, viability)

if (identical(colnames(sciPlex_A549_cds) , viability_column$cell)) {
  colData(sciPlex_A549_cds)$viability = viability_column$viability
} else{
  message("Make sure that rownames match before joining")
}


pData(sciPlex_A549_cds)$viability[is.na(pData(sciPlex_A549_cds)$viability)] <-
  1

sciPlex_A549_cds <- detect_genes(sciPlex_A549_cds, min_expr = 0.1)

sciPlex_A549_cds <- estimate_size_factors(sciPlex_A549_cds)

sciPlex_A549_cds <- estimate_cell_cycle(
  sciPlex_A549_cds,
  g1s_markers = cc.genes$s.genes,
  g2m_markers = cc.genes$g2m.genes
)

genes_expressed_per_drug =
  pData(sciPlex_A549_cds) %>%
  as.data.frame() %>%
  group_by(treatment, time_point) %>%
  nest() %>%
  mutate(fraction_genes = purrr::map(
    data,
    .f = function(pdata_subset, cds) {
      cds_subset = cds[, as.character(pdata_subset$cell)]
      cds_subset = detect_genes(cds_subset)
      tibble(
        id = fData(cds_subset)$id,
        gene_short_name = fData(cds_subset)$gene_short_name,
        num_cells_expressed = fData(cds_subset)$num_cells_expressed,
        fraction_cells_expressed = fData(cds_subset)$num_cells_expressed / ncol(cds_subset)
      )
    },
    sciPlex_A549_cds
  ))

genes_expressed_per_drug =
  genes_expressed_per_drug %>%
  unnest(fraction_genes)

expressed_genes = genes_expressed_per_drug %>%
  filter(fraction_cells_expressed > 0.05) %>%
  ungroup() %>%
  dplyr::select(id) %>%
  distinct()

expressed_genes = as.character(expressed_genes$id)



# After regresing out the effects of number of UMIs and replicate, perform PCA and return 25 PCs
sciPlex_A549_cds_umi_replicate <- preprocess_cds(
  sciPlex_A549_cds,
  method = 'PCA',
  num_dim = 25,
  norm_method = 'log',
  use_genes = expressed_genes,
  residual_model_formula_str = "~log(n.umi) + replicate",
  verbose = T,
  cores = 20
)

# Perform UMAP dimensionality reduction
sciPlex_A549_cds_umi_replicate <-
  reduce_dimension(
    sciPlex_A549_cds_umi_replicate,
    max_components = 2,
    reduction_method = 'UMAP',
    umap.metric = 'cosine',
    umap.n_neighbors = 10,
    umap.min_dist = 0.1,
    verbose = TRUE,
    cores = 20
  )

sciPlex_A549_cds_umi_replicate = append_umap_coordinates(sciPlex_A549_cds_umi_replicate)

ggplot(
  colData(sciPlex_A549_cds_umi_replicate) %>% as.data.frame(),
  aes(
    x = umap1,
    y = umap2,
    color = as.factor(time_point)
  )
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  scale_color_manual("Exposure\nTime", values = c("24" = "brown4", "72" =
                                                    "navy")) +
  ylim(-5, 7) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_umi_replicate_48_drug_subset.png",
    width = 2,
    height = 1,
    dpi = 900
  )



ggplot(
  colData(sciPlex_A549_cds_umi_replicate) %>% as.data.frame(),
  aes(x = umap1, y = umap2, color = g1s_score)
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  scale_color_viridis_c("G1/S \nScore") +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  ylim(-5, 7) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_umi_replicate_g1s_score.png",
    width = 2,
    height = 1,
    dpi = 900
  )


ggplot(
  colData(sciPlex_A549_cds_umi_replicate) %>% as.data.frame(),
  aes(x = umap1, y = umap2, color = g2m_score)
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  scale_color_viridis_c("G2/M \nScore") +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  ylim(-5, 7) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_umi_replicate_g2m_score.png",
    width = 2,
    height = 1,
    dpi = 900
  )


ggplot(
  colData(sciPlex_A549_cds_umi_replicate) %>% as.data.frame(),
  aes(x = umap1, y = umap2, color = g1s_score)
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  scale_color_viridis_c("G1/S \nScore") +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  ylim(-5, 7) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_umi_replicate_g1s_score.png",
    width = 2,
    height = 1,
    dpi = 900
  )


ggplot(
  colData(sciPlex_A549_cds_umi_replicate) %>% as.data.frame(),
  aes(x = umap1, y = umap2, color = proliferation_index)
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  scale_color_viridis_c("Proliferation \nIndex") +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  ylim(-5, 7) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_umi_replicate_PI_score.png",
    width = 2,
    height = 1,
    dpi = 900
  )



ggplot(
  colData(sciPlex_A549_cds_umi_replicate) %>% as.data.frame() %>% filter(vehicle),
  aes(x = umap1, y = umap2, color = proliferation_index)
) +
  geom_point(size = 0.05, stroke = 0, aes(color = vehicle)) +
  facet_wrap( ~ as.factor(time_point)) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  scale_color_manual("Vehicle",
                     values = c("TRUE" = "red"),
                     labels = c("Vehicle")) +
  ylim(-5, 7) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_umi_replicate_vehicle.png",
    width = 2,
    height = 1,
    dpi = 900
  )




# After regresing out the effects of number of UMIs and replicate, perform PCA and return 25 PCs
sciPlex_A549_cds_u_r_p_v <- preprocess_cds(
  sciPlex_A549_cds,
  method = 'PCA',
  num_dim = 25,
  norm_method = 'log',
  use_genes = expressed_genes,
  residual_model_formula_str = "~log(n.umi) + replicate + viability + proliferation_index",
  verbose = T,
  cores = 20
)

# Perform UMAP dimensionality reduction
sciPlex_A549_cds_u_r_p_v <-
  reduce_dimension(
    sciPlex_A549_cds_u_r_p_v,
    max_components = 2,
    reduction_method = 'UMAP',
    umap.metric = 'cosine',
    umap.n_neighbors = 10,
    umap.min_dist = 0.1,
    verbose = TRUE,
    cores = 20
  )

sciPlex_A549_cds_u_r_p_v = append_umap_coordinates(sciPlex_A549_cds_u_r_p_v)



ggplot(
  colData(sciPlex_A549_cds_u_r_p_v) %>% as.data.frame(),
  aes(
    x = umap1,
    y = umap2,
    color = as.factor(time_point)
  )
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  scale_color_manual("Exposure\nTime", values = c("24" = "brown4", "72" =
                                                    "navy")) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  ylim(-5, 8) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_urmpv_24_48.png",
    width = 2,
    height = 1,
    dpi = 900
  )


ggplot(
  colData(sciPlex_A549_cds_u_r_p_v) %>% as.data.frame(),
  aes(x = umap1, y = umap2, color = g2m_score)
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  scale_color_viridis_c("G2/M \nScore") +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  ylim(-5, 8) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_urmpv_g2m_score.png",
    width = 2,
    height = 1,
    dpi = 900
  )


ggplot(
  colData(sciPlex_A549_cds_u_r_p_v) %>% as.data.frame(),
  aes(x = umap1, y = umap2, color = g1s_score)
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  scale_color_viridis_c("G1/S \nScore") +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  ylim(-5, 8) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_urmpv_g1s_score.png",
    width = 2,
    height = 1,
    dpi = 900
  )


ggplot(
  colData(sciPlex_A549_cds_u_r_p_v) %>% as.data.frame(),
  aes(x = umap1, y = umap2, color = proliferation_index)
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  scale_color_viridis_c("Proliferation \nIndex") +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  ylim(-5, 8) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_urmpv_PI_score.png",
    width = 2,
    height = 1,
    dpi = 900
  )



ggplot(
  colData(sciPlex_A549_cds_u_r_p_v) %>% as.data.frame() %>% filter(vehicle),
  aes(x = umap1, y = umap2, color = proliferation_index)
) +
  geom_point(size = 0.05, stroke = 0, aes(color = vehicle)) +
  facet_wrap( ~ as.factor(time_point)) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  scale_color_manual("Vehicle",
                     values = c("TRUE" = "red"),
                     labels = c("Vehicle")) +
  ylim(-5, 8) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle:::monocle_theme_opts()
ggsave(
  "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_urmpv_vehicle.png",
  width = 2,
  height = 1,
  dpi = 900
)




cells_list <- list()
time_point <-
  unique(colData(sciPlex_A549_cds_u_r_p_v)[, "time_point"])

for (t in time_point) {
  cells_list[[as.character(t)]] <-
    row.names(pData(sciPlex_A549_cds_u_r_p_v)[pData(sciPlex_A549_cds_u_r_p_v)[, "time_point"] == t, ])
}

PCA_cells_1 = sciPlex_A549_cds_u_r_p_v@reducedDims[["PCA"]][cells_list[[1]], ]
PCA_cells_2 = sciPlex_A549_cds_u_r_p_v@reducedDims[["PCA"]][cells_list[[2]], ]



corrected_PCA <- scran::mnnCorrect(t(PCA_cells_1),
                                   t(PCA_cells_2))

corrected_PCA <- t(do.call("cbind", corrected_PCA[["corrected"]]))

row.names(corrected_PCA) <- c(row.names(PCA_cells_1),
                              row.names(PCA_cells_2))

sciPlex_A549_cds_urpv_MNN = sciPlex_A549_cds_u_r_p_v
corrected_PCA <- corrected_PCA[colnames(sciPlex_A549_cds_urpv_MNN), ]
sciPlex_A549_cds_urpv_MNN@reducedDims[["PCA"]] = as.matrix(corrected_PCA)




# Perform UMAP dimensionality reduction
sciPlex_A549_cds_urpv_MNN <-
  reduce_dimension(
    sciPlex_A549_cds_urpv_MNN,
    max_components = 2,
    reduction_method = 'UMAP',
    umap.metric = 'cosine',
    umap.n_neighbors = 10,
    umap.min_dist = 0.1,
    verbose = TRUE,
    cores = 20
  )

sciPlex_A549_cds_urpv_MNN = append_umap_coordinates(sciPlex_A549_cds_urpv_MNN)



ggplot(
  colData(sciPlex_A549_cds_urpv_MNN) %>% as.data.frame(),
  aes(
    x = umap1,
    y = umap2,
    color = as.factor(time_point)
  )
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  scale_color_manual("Exposure\nTime", values = c("24" = "brown4", "72" =
                                                    "navy")) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  ylim(-5, 7) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_urmpv_mnn_24_48.png",
    width = 2,
    height = 1,
    dpi = 900
  )

ggplot(
  colData(sciPlex_A549_cds_urpv_MNN) %>% as.data.frame(),
  aes(x = umap1, y = umap2, color = g1s_score)
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  scale_color_viridis_c("G1/S \nScore") +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  ylim(-5, 7) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_urmpv_mnn_g1s_score.png",
    width = 2,
    height = 1,
    dpi = 900
  )


ggplot(
  colData(sciPlex_A549_cds_urpv_MNN) %>% as.data.frame(),
  aes(x = umap1, y = umap2, color = g2m_score)
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  scale_color_viridis_c("G2/M \nScore") +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  ylim(-5, 7) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_urmpv_mnn_g2m_score.png",
    width = 2,
    height = 1,
    dpi = 900
  )

ggplot(
  colData(sciPlex_A549_cds_urpv_MNN) %>% as.data.frame(),
  aes(x = umap1, y = umap2, color = g1s_score)
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  scale_color_viridis_c("G1/S \nScore") +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  ylim(-5, 7) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_urmpv_mnn_g1s_score.png",
    width = 2,
    height = 1,
    dpi = 900
  )


ggplot(
  colData(sciPlex_A549_cds_urpv_MNN) %>% as.data.frame(),
  aes(x = umap1, y = umap2, color = proliferation_index)
) +
  geom_point(size = 0.05, stroke = 0) +
  facet_wrap( ~ as.factor(time_point)) +
  scale_color_viridis_c("Proliferation \nIndex") +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  ylim(-5, 7) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_urmpv_mnn_PI_score.png",
    width = 2,
    height = 1,
    dpi = 900
  )




ggplot(
  colData(sciPlex_A549_cds_urpv_MNN) %>% as.data.frame() %>% filter(vehicle),
  aes(x = umap1, y = umap2, color = proliferation_index)
) +
  geom_point(size = 0.05, stroke = 0, aes(color = vehicle)) +
  facet_wrap( ~ as.factor(time_point)) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.3, "line"),
    legend.key.height = unit(0.4, "line"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  ) +
  scale_color_manual("Vehicle",
                     values = c("TRUE" = "red"),
                     labels = c("Vehicle")) +
  ylim(-5, 7) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure24_sciPlex_A549_24-sciPlex_A549_24-72hr_UMAP_expressed_genes_by_time_point_Residual_model_urmpv_mnn_vehicle.png",
    width = 2,
    height = 1,
    dpi = 900
  )



ggplot(
  colData(sciPlex_A549_cds_urpv_MNN) %>% as.data.frame(),
  aes(x = umap1, y = umap2, color = pathway_level_1)
) +
  geom_point(size = 0.1, stroke = 0) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line")
  ) +
  scale_color_manual(
    "Pathway",
    values = c(
      "Antioxidant" = "aquamarine3",
      "Apoptotic regulation" = "lemonchiffon4",
      "Cell cycle regulation" = "deepskyblue2",
      "DNA damage & DNA repair" = "slategray4",
      "Epigenetic regulation" = "navy",
      "Focal adhesion signaling" =
        "brown4",
      "HIF signaling" = "darkgreen",
      "JAK/STAT signaling" = "orangered2",
      "Metabolic regulation" = "gold2",
      "Neuronal signaling" = "darkolivegreen",
      "Nuclear receptor signaling" = "chartreuse3",
      "PKC signaling" = "plum4",
      "Protein folding & Protein degradation" = "darkorchid2",
      "TGF/BMP signaling" = "darkcyan",
      "Tyrosine kinase signaling" = "firebrick1",
      "Other" = "darkgoldenrod4",
      "Vehicle" = "grey80"
    )
  ) +
  scale_alpha(aes(alpha = dose), range = c(0.01, 1)) +
  xlab("Component 1") +
  ylab("Component 2") +
  ylim(-5, 7) +
  facet_wrap( ~ as.factor(time_point)) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1.5)))) +
  monocle:::monocle_theme_opts() #+
ggsave(
  "supplemental_figure24_sciPlex_A549_24hr_UMAP_by_pathway_level_1_Residual_model_umi_replicate_vianbility_proliferation_index_time_point_aligned.png",
  width = 3.25,
  height = 2,
  unit = "in",
  dpi = 600
)

colData(sciPlex_A549_cds_urpv_MNN) %>%
  as.data.frame() %>%
  group_by(time_point) %>%
  add_tally() %>%
  ungroup() %>%
  group_by(time_point, pathway_level_1) %>%
  dplyr::mutate(fraction = n() / n) %>%
  dplyr::select(fraction, pathway_level_1) %>%
  distinct() %>%
  ggplot() +
  geom_bar(
    aes(
      x = reorder(pathway_level_1, fraction),
      y = fraction,
      fill = pathway_level_1
    ),
    stat = "identity",
    color = "black",
    size = .25
  ) +
  coord_flip() +
  facet_wrap( ~ as.factor(time_point)) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    strip.text = element_blank(),
    axis.ticks.x.top = element_line(),
    axis.line.x.top = element_line(),
    axis.text.x = element_text(angle = 90)
  ) +
  monocle:::monocle_theme_opts() +
  scale_y_continuous(position = "right", expand = c(0, 0)) +
  xlab("") +
  ylab("") +
  scale_fill_manual(
    "Pathway",
    values = c(
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
      "Other" = "darkgoldenrod4",
      "Vehicle" = "grey80"
    )
  ) +
  ggsave(
    "supplemental_figure24_Compounds_by_pathway_level_1.pdf",
    width = 2.5,
    height = 1.5,
    unit = "in",
    dpi = 900
  )


colData(sciPlex_A549_cds_urpv_MNN) %>%
  as.data.frame() %>%
  filter(pathway_level_1 == "Epigenetic regulation") %>%
  group_by(time_point) %>%
  add_tally() %>%
  ungroup() %>%
  group_by(time_point, pathway_level_2) %>%
  dplyr::mutate(fraction = n() / n) %>%
  dplyr::select(fraction, pathway_level_2, time_point) %>%
  distinct() %>%
  ggplot() +
  geom_bar(
    aes(
      x = reorder(pathway_level_2, fraction),
      y = fraction,
      fill = pathway_level_2
    ),
    stat = "identity",
    color = "black",
    size = .25
  ) +
  coord_flip() +
  facet_wrap( ~ as.factor(time_point)) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    strip.text = element_blank(),
    axis.ticks.x.top = element_line(),
    axis.line.x.top = element_line(),
    axis.text.x = element_text(angle = 90)
  ) +
  monocle:::monocle_theme_opts() +
  scale_y_continuous(position = "right", expand = c(0, 0)) +
  xlab("") +
  ylab("") +
  scale_fill_manual(
    "Pathway",
    values  = c(
      "JAK kinase activity" = "aquamarine3",
      "RTK activity" = "lemonchiffon4",
      "Histone deacetylation" = "deepskyblue2",
      "Mitochondria-mediated apoptosis" = "slategray4",
      "Histone methylation" = "firebrick1",
      "HSP90 activity" = "brown4",
      "MAPK activity" = "darkgreen",
      "DNA methylation" = "orangered3",
      "Cell cycle regulation" = "gold4",
      "Bromodomain" = "darkolivegreen",
      "CDK activity" = "chartreuse3",
      "PIM1 activity" = "plum4",
      "Nuclear receptor activity" = "darkorchid2",
      "PKC activitiy" = "darkcyan",
      "GPCR activity" = "navy",
      "ADP-rybosilation" = "darkgoldenrod4",
      "Nucleotide analog" = "green4",
      "Toposiomerase activity" = "cyan3",
      "Histamine receptor activity" = "bisque3",
      "PI3K-AKT-MTOR activity" = "darkslategray2",
      "Isocitrate dehydrogenase activity" = "palevioletred2"
    )
  ) +
  ggsave(
    "supplemental_figure24_Epigenetic_compounds_by_pathway_level_2.pdf",
    width = 2.5,
    height = 1.5,
    unit = "in",
    dpi = 900
  )


colData(sciPlex_A549_cds_urpv_MNN) %>%
  as.data.frame() %>%
  filter(pathway_level_2 == "Histone deacetylation") %>%
  group_by(time_point) %>%
  add_tally() %>%
  ungroup() %>%
  group_by(time_point, product_name) %>%
  dplyr::mutate(fraction = n() / n) %>%
  dplyr::select(fraction, product_name, time_point) %>%
  distinct() %>%
  ggplot() +
  geom_bar(
    aes(x = reorder(product_name, fraction), y = fraction),
    stat = "identity",
    color = "black",
    size = .25,
    fill = "grey80"
  ) +
  coord_flip() +
  facet_wrap( ~ as.factor(time_point)) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    strip.text = element_blank(),
    axis.ticks.x.top = element_line(),
    axis.line.x.top = element_line(),
    axis.text.x = element_text(angle = 90)
  ) +
  monocle:::monocle_theme_opts() +
  scale_y_continuous(position = "right", expand = c(0, 0)) +
  xlab("") +
  ylab("") +
  ggsave(
    "supplemental_figure24_HDAC_compounds_by_product.pdf",
    width = 2.5,
    height = 1.5,
    unit = "in",
    dpi = 900
  )

sciPlex_A549_cds_urpv_MNN %>%
  colData() %>%
  as.data.frame() %>%
  filter(product_name %in% c("Vehicle", "Belinostat (PXD101)")) %>%
  ggplot(aes(x = umap1, y = umap2, color = as.factor(dose))) +
  geom_point(size = 0.25, stroke = 0) +
  facet_wrap( ~ time_point) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    strip.text = element_blank()
  ) +
  scale_color_manual(
    "Log(Dose Belinostat [nM])",
    values = c(
      "0" = "gray",
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure25_supplemental_figure25_Aligned_24-72hr_A549_UMAP_proliferation_index_by_belinostat_dose.png",
    width = 3.5,
    height = 2
  )


sciPlex_A549_cds_urpv_MNN %>%
  colData() %>%
  as.data.frame() %>%
  filter(product_name %in% c("Vehicle", "Abexinostat (PCI-24781)")) %>%
  ggplot(aes(x = umap1, y = umap2, color = as.factor(dose))) +
  geom_point(size = 0.25, stroke = 0) +
  facet_wrap( ~ time_point) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    strip.text = element_blank()
  ) +
  scale_color_manual(
    "Log(Abexinostat [nM])",
    values = c(
      "0" = "gray",
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure25_Aligned_24-72hr_A549_UMAP_proliferation_index_by_Abexinostat_dose.png",
    width = 3.5,
    height = 2
  )

sciPlex_A549_cds_urpv_MNN %>%
  colData() %>%
  as.data.frame() %>%
  filter(product_name %in% c("Vehicle", "(+)-JQ1")) %>%
  ggplot(aes(x = umap1, y = umap2, color = as.factor(dose))) +
  geom_point(size = 0.25, stroke = 0) +
  facet_wrap( ~ time_point) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    strip.text = element_blank()
  ) +
  scale_color_manual(
    "Log(Dose JQ1 [nM])",
    values = c(
      "0" = "gray",
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure25_Aligned_24-72hr_A549_UMAP_proliferation_index_by_JQ1_dose.png",
    width = 3.5,
    height = 2
  )


sciPlex_A549_cds_urpv_MNN %>%
  colData() %>%
  as.data.frame() %>%
  filter(product_name %in% c("Vehicle", "SRT2104 (GSK2245840)")) %>%
  ggplot(aes(x = umap1, y = umap2, color = as.factor(dose))) +
  geom_point(size = 0.25, stroke = 0) +
  facet_wrap( ~ time_point) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    strip.text = element_blank()
  ) +
  scale_color_manual(
    "Log(Dose SRT2104 [nM])",
    values = c(
      "0" = "gray",
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure25_Aligned_24-72hr_A549_UMAP_proliferation_index_by_SRT2104_dose.png",
    width = 3.5,
    height = 2
  )
