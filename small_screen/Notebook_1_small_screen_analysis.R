
# Input path to github directory
path_to_github = "~/sci-Plex/"
path_to_monocle3 = paste0(path_to_github,
                          "monocle3_d4a9a35/monocle3/",
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
  library(devtools)
  load_all(path_to_monocle3)
})

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)

setwd(paste0(path_to_github,
             "small_screen",
             sep = ""))

# Set directory for sciPlex github 
bin_directory = paste0(path_to_github,
                       "bin",
                       sep = "")

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
  cc.genes <- readRDS(paste0(bin_directory,
                             "cc.genes.RDS",
                             sep = ""))
  source(paste0(bin_directory,
                "viability.R",
                sep = ""))
  source(paste0(bin_directory,
                "dispersions_functions.R",
                sep = ""))
  
})


# Read in CDS file from experiment
# CDS can be found at GEO
sciPlex_cds <- readRDS("sciPlex2_cds.RDS")

sciPlex_cds <-
  sciPlex_cds[, !is.na(pData(sciPlex_cds)$top_to_second_best_ratio)]

dim(sciPlex_cds)

ggplot(pData(sciPlex_cds) %>% as.data.frame(),
       aes(x = log10(top_to_second_best_ratio))) +
  geom_histogram() +
  monocle3::monocle_theme_opts()

ggplot(pData(sciPlex_cds)[pData(sciPlex_cds)$top_to_second_best_ratio > 10 &
                            pData(sciPlex_cds)$qval < 0.01 &
                            pData(sciPlex_cds)$hash_umis > 30, ]  %>% as.data.frame(), aes(x = log10(top_to_second_best_ratio))) +
  geom_histogram() +
  monocle3::monocle_theme_opts()

# Subset cells that have more than a 10 fold enrichment in hash oligos, and over 30 hash oligo UMIs
sciPlex_cds <-
  sciPlex_cds[, pData(sciPlex_cds)$top_to_second_best_ratio > 10 &
                pData(sciPlex_cds)$qval < 0.01 &
                pData(sciPlex_cds)$hash_umis > 30]

colData(sciPlex_cds)$cell_type = "A549"

meta_data = 
  stringr::str_split_fixed(colData(sciPlex_cds)$top_oligo,
                                     pattern = "_",
                                     n = 3)

colData(sciPlex_cds)$cell =
  colData(sciPlex_cds)$Cell

colData(sciPlex_cds)$dose =
  meta_data[, 2] %>%
  as.numeric()

colData(sciPlex_cds)$treatment =
  meta_data[, 1] %>%
  as.character()

colData(sciPlex_cds)$culture_well =
  meta_data[, 3] %>%
  stringr::str_sub(start = 2,
                   end = 4) %>%
  as.character()

colData(sciPlex_cds)$well_oligo =
  meta_data[, 3] %>%
  as.character()


colData(sciPlex_cds)$culture_plate =
  meta_data[, 3] %>%
  stringr::str_sub(start = 1,
                   end = 1) %>%
  as.character()

colData(sciPlex_cds)$vehicle =
  colData(sciPlex_cds)$dose  == 0

colData(sciPlex_cds)$dose_character =
  factor(colData(sciPlex_cds)$dose,
         levels = c("0", "0.1", "0.5", "1", "5", "10", "50", "100"))


colData(sciPlex_cds)$new_treatment_label =
  sapply(colData(sciPlex_cds)$treatment, function(x) {
    if (grepl("BMS", x))
      return("BMS345541")
    if (grepl("Dex", x))
      return("Dex")
    if (grepl("Nutlin", x))
      return("Nutlin3A")
    if (grepl("SAHA", x))
      return("SAHA")
  })


colData(sciPlex_cds)$solvent =
  sapply(pData(sciPlex_cds)$treatment, function(x) {
    if (grepl("Dex", x))
      return("Ethanol")
    else
      return("DMSO")
  })


colData(sciPlex_cds)$replicate <-
  unlist(sapply(colData(sciPlex_cds)$culture_well, function(x) {
    if (stringr::str_sub(x, -2) %in% c("01", "04", "07", "10"))
      return("rep1")
    if (stringr::str_sub(x, -2) %in% c("02", "05", "08", "11"))
      return("rep2")
    if (stringr::str_sub(x, -2) %in% c("03", "06", "09", "12"))
      return("rep3")
    else
      return(NA)
  }))

# Calculate dose response relationships for each drug treatment
sciPlex_counts = count_cells(sciPlex_cds)
sciPlex_models = fit_dose_response_models(sciPlex_counts)
sciPlex_model_stats = extract_dose_model_stats(sciPlex_counts, sciPlex_models)

# Calculate the proliferation index for every cell
sciPlex_cds = estimate_viability(sciPlex_cds, sciPlex_models)

# Pre-process cds
sciPlex_cds = detect_genes(sciPlex_cds, min_expr = 0.1)
sciPlex_cds = estimate_size_factors(sciPlex_cds)

expressed_genes = row.names(fData(sciPlex_cds)[Matrix::rowSums(as.matrix(exprs(sciPlex_cds)) > 0) > 50 , ])
length(expressed_genes)

# Identify over-dispersed genes to be used as feature selection
dispersion_test = estimateDispersionsForCellDataSet(sciPlex_cds,
                                                    min_cells_detected = 100)

disp_table = dispersionTable(dispersion_test)

disp_table = disp_table %>%
  mutate(excess_disp = (dispersion_empirical - dispersion_fit) / dispersion_fit) %>%
  dplyr::arrange(plyr::desc(excess_disp))

unsup_clustering_genes = as.character(head(disp_table, 1000)$gene_id)

ordering_genes = unsup_clustering_genes


# Perform dimensionality redction using UMAP
sciPlex_cds <- monocle3::preprocess_cds(
  sciPlex_cds,
  method = 'PCA',
  num_dim = 35,
  norm_method = 'log',
  residualModelFormulaStr = "~log(n.umi)",
  use_genes = ordering_genes,
  verbose = T
)



sciPlex_cds <- reduce_dimension(
  sciPlex_cds,
  max_components = 2,
  reduction_method = 'UMAP',
  umap.metric = "cosine",
  umap.n_neighbors = 50,
  umap.min_dist = 0.5,
  umap.fast_sgd = FALSE,
  cores = 1,
  verbose = T
)

# Cluster cells based on PCA
sciPlex_cds = cluster_cells(sciPlex_cds,
                            reduction_method = "PCA")

pData(sciPlex_cds)$Cluster = clusters(sciPlex_cds,
                                      reduction_method = "PCA")

# Copy UMAP coordinates to colData
colData(sciPlex_cds)$UMAP_1 = sciPlex_cds@reducedDims[["UMAP"]][, 1]
colData(sciPlex_cds)$UMAP_2 = sciPlex_cds@reducedDims[["UMAP"]][, 2]

# Plots
pdf(
  "supplemental_figure2_sciPlex2_UMAP_by_PCA_cluster.pdf",
  width = 2.5,
  height = 2
)
ggplot(pData(sciPlex_cds) %>% as.data.frame(),
       aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
  geom_point(size = 0.3, stroke = 0) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line")
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle3:::monocle_theme_opts()
dev.off()

pdf("figure2_sciPlex2_UMAP_by_treatment.pdf",
    width = 2,
    height = 1.25)
ggplot(pData(sciPlex_cds) %>% as.data.frame(),
       aes(x = UMAP_1, y = UMAP_2, color = treatment)) +
  geom_point(size = 0.2, stroke = 0) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.25, "line")
  ) +
  scale_color_manual(
    "Treatment",
    labels = c(
      "BMS" = "BMS345541",
      "Dex" = "Dex",
      "Nutlin" = "Nutlin3A",
      "SAHA" = "SAHA"
    ),
    values = c(
      "BMS" = "dimgrey",
      "Dex" = "deepskyblue3",
      "Nutlin" = "firebrick3",
      "SAHA" = "springgreen4"
    )
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1.5)))) +
  monocle3:::monocle_theme_opts()
dev.off()

vehicle_subset_pData <-
  pData(sciPlex_cds)[pData(sciPlex_cds)$vehicle == TRUE, c("UMAP_1", "UMAP_2")]

pdf("figure2_sciPlex2_UMAP_by_log_dose.pdf",
    width = 2,
    height = 1.75)
ggplot(
  pData(sciPlex_cds) %>% as.data.frame(),
  aes(x = UMAP_1, y = UMAP_2, color = dose_character)
) +
  geom_point(
    data = vehicle_subset_pData  %>% as.data.frame(),
    color = "gray",
    size = 0.05,
    stroke = 0
  ) +
  geom_point(size = 0.2, stroke = 0) +
  facet_wrap( ~ new_treatment_label, ncol = 2) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.25, "line")
  ) +
  scale_color_manual(
    "Dose",
    values = c(
      "0" = "gray",
      "0.1" = "#1D1147FF",
      "0.5" = "#51127CFF",
      "1" = "#822681FF",
      "5" = "#B63679FF",
      "10" = "#E65164FF",
      "50" = "#FB8861FF",
      "100" = "#FEC287FF"
    )
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1.5)))) +
  monocle3:::monocle_theme_opts()
dev.off()

pdf(
  "supplemental_figure2_sciPlex2_UMAP_by_vehicle.pdf",
  width = 2.5,
  height = 2
)
ggplot(pData(sciPlex_cds) %>% as.data.frame(),
       aes(x = UMAP_1, y = UMAP_2, color = vehicle)) +
  geom_point(size = 0.3, stroke = 0) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.25, "line")
  ) +
  scale_color_manual(
    "Treatment",
    labels = c("TRUE" = "Vehicle", "FALSE" = "Treated"),
    values = c("FALSE" = "gray", "TRUE" = "firebrick4")
  ) +
  geom_point(
    data = pData(sciPlex_cds)[pData(sciPlex_cds)$vehicle == "TRUE", ] %>% as.data.frame(),
    aes(x = UMAP_1, y = UMAP_2),
    color = "firebrick4",
    size = 0.05,
    stroke = 0
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle3:::monocle_theme_opts()
dev.off()

pdf("figure2_sciPlex2_dose_response_curves.pdf",
    width = 2,
    height = 2)
plot_dose_response_curves(sciPlex_counts,
                          sciPlex_models,
                          color_by = "factor(dose)",
                          ncol = 2) +
  xlab("Log(Dose [µM])") +
  ylab("Cell counts") +
  theme(
    text = element_text(size = 6),
    legend.position = "none",
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.5, "line")
  ) +
  scale_color_manual(
    "Dose",
    values = c(
      "0" = "gray",
      "0.1" = "#1D1147FF",
      "0.5" = "#51127CFF",
      "1" = "#822681FF",
      "5" = "#B63679FF",
      "10" = "#E65164FF",
      "50" = "#FB8861FF",
      "100" = "#FEC287FF"
    )
  ) +
  facet_wrap( ~ treatment, scales = "free") +
  ylim(0, NA)
dev.off()


# Compare cell counts obtained via sci-Plex to cell counts obtained via cell titer glo
sciPlex_counts <-
  sciPlex_counts %>% group_by(treatment) %>% mutate(vehicle_cells = mean(norm_cells[dose == 0]))

sciPlex_counts$vehicle_norm_cells <-
  sciPlex_counts$norm_cells / sciPlex_counts$vehicle_cells

norm_sciPlex_counts <- as.data.frame(
  sciPlex_counts %>%
    #filter(dose %in% c(10, 5, 1, 0.5, 0.1, 0)) %>%
    group_by(dose, treatment) %>%
    dplyr::summarise(
      mean_counts = mean(vehicle_norm_cells),
      stdev_counts = sd(vehicle_norm_cells)
    )
)

row.names(norm_sciPlex_counts) <-
  paste0(norm_sciPlex_counts$treatment, "_", norm_sciPlex_counts$dose)

ctg_values <- read.csv("sciPlex2_CTG_results.csv")


norm_ctg_values <-
  as.data.frame(
    ctg_values %>% group_by(dose, treatment) %>%
      dplyr::summarise(mean_CTG = mean(CTG), stdev_CTG = sd(CTG))
  )

row.names(norm_ctg_values) <-
  paste0(norm_ctg_values$treatment, "_", norm_ctg_values$dose)

dose_response_values_df <-
  merge(norm_sciPlex_counts,
        norm_ctg_values %>% dplyr::select(mean_CTG, stdev_CTG),
        by = "row.names")


pdf(
  "supplemental_figure2_sciPlex2_Hash_counts_vs_CTG_FigS2.pdf",
  width = 2,
  height = 1.75
)
ggplot(dose_response_values_df) +
  geom_point(
    aes(
      x = log(dose + 0.1),
      y = mean_counts,
      color = "Nuclear Hash counts"
    ),
    color = "dimgrey",
    size = 1,
    stroke = 0
  ) +
  geom_point(
    aes(
      x = log(dose + 0.1),
      y = mean_CTG,
      color = "Cell-titer Glo"
    ),
    color = "brown4",
    size = 1,
    stroke = 0
  ) +
  facet_wrap( ~ treatment) +
  ylab("Relative viability") +
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    axis.text.x = element_text(size = 4)
  )
dev.off()

cor.test(dose_response_values_df$mean_counts,
         dose_response_values_df$mean_CTG,
         method = "pearson")

ggplot(dose_response_values_df, aes(x = mean_counts, y = mean_CTG)) +
  geom_point(size = 1, stroke = 0) +
  geom_errorbar(aes(ymin = mean_CTG - stdev_CTG, ymax = mean_CTG + stdev_CTG), lwd = 0.1) +
  geom_errorbarh(aes(xmin = mean_counts - stdev_counts, xmax = mean_counts + stdev_counts),
                 lwd = 0.1) +
  #facet_wrap(~treatment, scales = "free") +
  geom_smooth(method = "lm") +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    axis.text.x = element_text(size = 4)
  ) +
  xlim(0, NA) +
  ylim(0, NA) +
  xlab("Inferred cell counts\n(Mean value relative to vehicle)") +
  ylab("Cell titer glo\n(Mean value relative to vehicle)") +
  monocle_theme_opts() +
  ggsave(
    "supplemental_figure2_Correlation_between_inferred_cell_counts_vs_CTG.png",
    width = 2,
    height = 1.75,
    dpi = 600
  )

# Plot markers of p53 and GR activity
apoptosis_markers <-
  sciPlex_cds[fData(sciPlex_cds)$gene_short_name %in% c("PMAIP1", "BBC3"),
              pData(sciPlex_cds)$treatment == "Nutlin" |
                pData(sciPlex_cds)$vehicle == TRUE]

pdf(
  "figure2_Apoptosis_markers_percent_positive_cells_nutlin_treated_cells.pdf",
  width = 2,
  height = 0.8
)
plot_percent_cells_positive(apoptosis_markers,
                            group_cells_by = "dose_character",
                            ncol = 2) +
  scale_fill_manual(
    "Dose (µM)",
    values = c(
      "0" = "gray",
      "0.1" = "#1D1147FF",
      "0.5" = "#51127CFF",
      "1" = "#822681FF",
      "5" = "#B63679FF",
      "10" = "#E65164FF",
      "50" = "#FB8861FF",
      "100" = "#FEC287FF"
    )
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 4)
  ) +
  xlab("Log(Dose [µM])") +
  ylab("Cells (percent)")
dev.off()

p53_nutlin_markers <-
  sciPlex_cds[fData(sciPlex_cds)$gene_short_name %in% c("CDKN1A", "TP53I3"),
              pData(sciPlex_cds)$treatment == "Nutlin" |
                pData(sciPlex_cds)$vehicle == TRUE]

pdf(
  "figure2_p53_marker_percent_positive_cells_nutlin_treated_cells.pdf",
  width = 2,
  height = 0.8
)
plot_percent_cells_positive(
  p53_nutlin_markers,
  group_cells_by = "dose_character",
  ncol = 2,
  plot_as_count = FALSE
) +
  scale_fill_manual(
    "Dose (µM)",
    values = c(
      "0" = "gray",
      "0.1" = "#1D1147FF",
      "0.5" = "#51127CFF",
      "1" = "#822681FF",
      "5" = "#B63679FF",
      "10" = "#E65164FF",
      "50" = "#FB8861FF",
      "100" = "#FEC287FF"
    )
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 4)
  ) +
  xlab("Log(Dose [µM])") +
  ylab("Cells (percent)")
dev.off()

Dex_markers <-
  sciPlex_cds[fData(sciPlex_cds)$gene_short_name %in% c("ANGPTL4", "GDF15"),
              pData(sciPlex_cds)$treatment == "Dex" |
                pData(sciPlex_cds)$vehicle == TRUE]

pdf(
  "figure2_GR_marker_percent_positive_cells_dex_treated_cells.pdf",
  width = 2,
  height = 0.8
)
plot_percent_cells_positive(
  Dex_markers,
  group_cells_by = "dose_character",
  ncol = 2,
  plot_as_count = FALSE
) +
  scale_fill_manual(
    "Dose (µM)",
    values = c(
      "0" = "gray",
      "0.1" = "#1D1147FF",
      "0.5" = "#51127CFF",
      "1" = "#822681FF",
      "5" = "#B63679FF",
      "10" = "#E65164FF",
      "50" = "#FB8861FF",
      "100" = "#FEC287FF"
    )
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 4)
  ) +
  xlab("Log(Dose [µM])") +
  ylab("Cells (percent)")
dev.off()
