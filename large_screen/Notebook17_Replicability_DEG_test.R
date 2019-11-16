# devtools::install_github(repo = 'cole-trapnell-lab/monocle3', ref = "d4a9a35")

# Input path to github directory
path_to_github = "~/sci-Plex/"

# Directory of sciPlex bin
bin_directory = paste0(path_to_github,
                       "bin",
                       sep = "")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(Matrix)
  library(devtools)
  library(monocle3)
})

# Load clean cds objects that include cells from rescue phenocopy sciPlex experiments and large screen
cds <- readRDS("sciPlex4.RDS")

# Create cds subsets where each subset is composed of cells exposed to
# 10 µM Pracinostat, Abexinostat and DMSO for sciPlex20 and the 3 level screen

A549_cds.list <- list()

A549_cds.list[["sciPlex20"]] <- list()
A549_cds.list[["Screen"]] <- list()


A549_cds.list[["Screen"]][["Abexinostat"]] <-
  cds[, (
    pData(cds)$new_hdac == "Abexinostat" &
      pData(cds)$hdac_dose == 10000 &
      pData(cds)$sample == "Plex3_Screen2" &
      pData(cds)$cell_type == "A549"
  ) |
    (
      pData(cds)$new_hdac == "Vehicle" &
        pData(cds)$sample == "Plex3_Screen2" &
        pData(cds)$cell_type == "A549"
    )]

A549_cds.list[["Screen"]][["Pracinostat"]] <-
  cds[, (
    pData(cds)$new_hdac == "Pracinostat" &
      pData(cds)$hdac_dose == 10000 &
      pData(cds)$sample == "Plex3_Screen2" &
      pData(cds)$cell_type == "A549"
  ) |
    (
      pData(cds)$new_hdac == "Vehicle" &
        pData(cds)$sample == "Plex3_Screen2" &
        pData(cds)$cell_type == "A549"
    )]

A549_cds.list[["sciPlex20"]][["Abexinostat"]] <-
  cds[, (
    pData(cds)$new_hdac == "Abexinostat" &
      pData(cds)$hdac_dose == 10000 &
      pData(cds)$treatment == "DMSO" &
      pData(cds)$sample == "sciPlex20" &
      pData(cds)$cell_type == "A549"
  ) |
    (
      pData(cds)$new_hdac == "Vehicle" &
        pData(cds)$treatment == "DMSO" &
        pData(cds)$sample == "sciPlex20" &
        pData(cds)$cell_type == "A549"
    )]

A549_cds.list[["sciPlex20"]][["Pracinostat"]] <-
  cds[, (
    pData(cds)$new_hdac == "Pracinostat" &
      pData(cds)$hdac_dose == 10000 &
      pData(cds)$treatment == "DMSO" &
      pData(cds)$sample == "sciPlex20" &
      pData(cds)$cell_type == "A549"
  ) |
    (
      pData(cds)$new_hdac == "Vehicle" &
        pData(cds)$treatment == "DMSO" &
        pData(cds)$sample == "sciPlex20" &
        pData(cds)$cell_type == "A549"
    )]

MCF7_cds.list <- list()

MCF7_cds.list[["sciPlex20"]] <- list()
MCF7_cds.list[["Screen"]] <- list()


MCF7_cds.list[["Screen"]][["Abexinostat"]] <-
  cds[, (
    pData(cds)$new_hdac == "Abexinostat" &
      pData(cds)$hdac_dose == 10000 &
      pData(cds)$sample == "Plex3_Screen2" &
      pData(cds)$cell_type == "MCF7"
  ) |
    (
      pData(cds)$new_hdac == "Vehicle" &
        pData(cds)$sample == "Plex3_Screen2" &
        pData(cds)$cell_type == "MCF7"
    )]

MCF7_cds.list[["Screen"]][["Pracinostat"]] <-
  cds[, (
    pData(cds)$new_hdac == "Pracinostat" &
      pData(cds)$hdac_dose == 10000 &
      pData(cds)$sample == "Plex3_Screen2" &
      pData(cds)$cell_type == "MCF7"
  ) |
    (
      pData(cds)$new_hdac == "Vehicle" &
        pData(cds)$sample == "Plex3_Screen2" &
        pData(cds)$cell_type == "MCF7"
    )]

MCF7_cds.list[["sciPlex20"]][["Abexinostat"]] <-
  cds[, (
    pData(cds)$new_hdac == "Abexinostat" &
      pData(cds)$hdac_dose == 10000 &
      pData(cds)$treatment == "DMSO" &
      pData(cds)$sample == "sciPlex20" &
      pData(cds)$cell_type == "MCF7"
  ) |
    (
      pData(cds)$new_hdac == "Vehicle" &
        pData(cds)$treatment == "DMSO" &
        pData(cds)$sample == "sciPlex20" &
        pData(cds)$cell_type == "MCF7"
    )]

MCF7_cds.list[["sciPlex20"]][["Pracinostat"]] <-
  cds[, (
    pData(cds)$new_hdac == "Pracinostat" &
      pData(cds)$hdac_dose == 10000 &
      pData(cds)$treatment == "DMSO" &
      pData(cds)$sample == "sciPlex20" &
      pData(cds)$cell_type == "MCF7"
  ) |
    (
      pData(cds)$new_hdac == "Vehicle" &
        pData(cds)$treatment == "DMSO" &
        pData(cds)$sample == "sciPlex20" &
        pData(cds)$cell_type == "MCF7"
    )]


# Perform differential gene expression analysis as a function of HDACi exposure for cells
# from our original screen and the rescue/phenocopy experiment

A549_screen_HDACi_diff_test_results <- list()
for (treatment in names(A549_cds.list[["Screen"]])) {
  A549_cds.list[["Screen"]][[treatment]] <-
    estimate_size_factors(A549_cds.list[["Screen"]][[treatment]])
  
  model_results <-
    fit_models(A549_cds.list[["Screen"]][[treatment]],
               model_formula_str = "~hdac_dose",
               cores = 40)
  
  A549_screen_HDACi_diff_test_results[[treatment]] <-
    coefficient_table(model_results)
  closeAllConnections()
}
saveRDS(A549_screen_HDACi_diff_test_results,
        "A549_screen_HDACi_diff_test_results.rds")


A549_sciplex20_HDACi_diff_test_results <- list()
for (treatment in names(A549_cds.list[["sciPlex20"]])) {
  A549_cds.list[["sciPlex20"]][[treatment]] <-
    estimate_size_factors(A549_cds.list[["sciPlex20"]][[treatment]])
  
  model_results <-
    fit_models(A549_cds.list[["sciPlex20"]][[treatment]],
               model_formula_str = "~hdac_dose",
               cores = 40)
  
  A549_sciplex20_HDACi_diff_test_results[[treatment]] <-
    coefficient_table(model_results)
  closeAllConnections()
}

saveRDS(
  A549_sciplex20_HDACi_diff_test_results,
  "A549_sciplex20_HDACi_diff_test_results.rds"
)

MCF7_screen_HDACi_diff_test_results <- list()
for (treatment in names(MCF7_cds.list[["Screen"]])) {
  MCF7_cds.list[["Screen"]][[treatment]] <-
    estimate_size_factors(MCF7_cds.list[["Screen"]][[treatment]])
  
  model_results <-
    fit_models(MCF7_cds.list[["Screen"]][[treatment]],
               model_formula_str = "~hdac_dose",
               cores = 40)
  
  MCF7_screen_HDACi_diff_test_results[[treatment]] <-
    coefficient_table(model_results)
  closeAllConnections()
}

saveRDS(MCF7_screen_HDACi_diff_test_results,
        "MCF7_screen_HDACi_diff_test_results.rds")

MCF7_sciplex20_HDACi_diff_test_results <- list()

for (treatment in names(MCF7_cds.list[["sciPlex20"]])) {
  MCF7_cds.list[["sciPlex20"]][[treatment]] <-
    estimate_size_factors(MCF7_cds.list[["sciPlex20"]][[treatment]])
  
  model_results <-
    fit_models(MCF7_cds.list[["sciPlex20"]][[treatment]],
               model_formula_str = "~hdac_dose",
               cores = 40)
  
  MCF7_sciplex20_HDACi_diff_test_results[[treatment]] <-
    coefficient_table(model_results)
  
  closeAllConnections()
}

saveRDS(MCF7_sciplex20_HDACi_diff_test_results,
        "MCF7_sciplex20_HDACi_diff_test_results.rds"
        )

# Make a list of DEGs as a function of HDACi epxosure across both experiments

# For Abexinostat
A549_Abexinostat_deg.list <- list()

A549_Abexinostat_deg.list[["Screen"]] <-
  (
    subset(
      A549_screen_HDACi_diff_test_results[["Abexinostat"]],
      q_value < 0.05 & term == "hdac_dose" & normalized_effect != 0
    )
  )$id

A549_Abexinostat_deg.list[["sciPlex20"]] <-
  (
    subset(
      A549_sciplex20_HDACi_diff_test_results[["Abexinostat"]],
      q_value < 0.05 & term == "hdac_dose" & normalized_effect != 0
    )
  )$id

A549_Abexinostat_degs <-
  unique(union(A549_Abexinostat_deg.list[["Screen"]], A549_Abexinostat_deg.list[["sciPlex20"]]))

# For Pracinostat
A549_Pracinostat_deg.list <- list()

A549_Pracinostat_deg.list[["Screen"]] <-
  (
    subset(
      A549_screen_HDACi_diff_test_results[["Pracinostat"]],
      q_value < 0.05 & term == "hdac_dose" & normalized_effect != 0
    )
  )$id

A549_Pracinostat_deg.list[["sciPlex20"]] <-
  (
    subset(
      A549_sciplex20_HDACi_diff_test_results[["Pracinostat"]],
      q_value < 0.05 & term == "hdac_dose" & normalized_effect != 0
    )
  )$id

A549_Pracinostat_degs <-
  unique(union(A549_Pracinostat_deg.list[["Screen"]], A549_Pracinostat_deg.list[["sciPlex20"]]))

# For Abexinostat
MCF7_Abexinostat_deg.list <- list()

MCF7_Abexinostat_deg.list[["Screen"]] <-
  (
    subset(
      MCF7_screen_HDACi_diff_test_results[["Abexinostat"]],
      q_value < 0.05 & term == "hdac_dose" & normalized_effect != 0
    )
  )$id

MCF7_Abexinostat_deg.list[["sciPlex20"]] <-
  (
    subset(
      MCF7_sciplex20_HDACi_diff_test_results[["Abexinostat"]],
      q_value < 0.05 & term == "hdac_dose" & normalized_effect != 0
    )
  )$id

MCF7_Abexinostat_degs <-
  unique(union(MCF7_Abexinostat_deg.list[["Screen"]], MCF7_Abexinostat_deg.list[["sciPlex20"]]))

# For Pracinostat
MCF7_Pracinostat_deg.list <- list()

MCF7_Pracinostat_deg.list[["Screen"]] <-
  (
    subset(
      MCF7_screen_HDACi_diff_test_results[["Pracinostat"]],
      q_value < 0.05 & term == "hdac_dose" & normalized_effect != 0
    )
  )$id

MCF7_Pracinostat_deg.list[["sciPlex20"]] <-
  (
    subset(
      MCF7_sciplex20_HDACi_diff_test_results[["Pracinostat"]],
      q_value < 0.05 & term == "hdac_dose" & normalized_effect != 0
    )
  )$id

MCF7_Pracinostat_degs <-
  unique(union(MCF7_Pracinostat_deg.list[["Screen"]], MCF7_Pracinostat_deg.list[["sciPlex20"]]))


# Create data frames of effect sizes and significance values for HDACi exposed cells
# across both experiments

# For Abexinostat
A549_screen_abexinostat_effect_sizes <-
  A549_screen_HDACi_diff_test_results[["Abexinostat"]] %>%
  dplyr::select(id, gene_short_name, q_value, normalized_effect)

colnames(A549_screen_abexinostat_effect_sizes) <-
  c("id",
    "gene_short_name",
    "q_value_screen",
    "normalized_effect_screen")

A549_sciplex20_abexinostat_effect_sizes <-
  A549_sciplex20_HDACi_diff_test_results[["Abexinostat"]] %>%
  dplyr::select(id, gene_short_name, q_value, normalized_effect)

colnames(A549_sciplex20_abexinostat_effect_sizes) <-
  c("id",
    "gene_short_name",
    "q_value_sciplex20",
    "normalized_effect_sciplex20")

# Ensure data frames are in the same order although left_join should merge safely
identical(
  A549_screen_abexinostat_effect_sizes$id,
  A549_sciplex20_abexinostat_effect_sizes$id
)

# Generate a joint data frame
A549_abexinostat_effect_sizes <-
  left_join(
    A549_screen_abexinostat_effect_sizes,
    A549_sciplex20_abexinostat_effect_sizes[, c("id", "q_value_sciplex20", "normalized_effect_sciplex20")],
    by = "id"
  )

A549_abexinostat_effect_sizes_subset <-
  A549_abexinostat_effect_sizes %>% filter(id %in% A549_Abexinostat_degs)
A549_abexinostat_effect_sizes_subset <-
  A549_abexinostat_effect_sizes_subset %>% filter(normalized_effect_screen != 0 &
                                                    normalized_effect_sciplex20 != 0)

cor.test(
  A549_abexinostat_effect_sizes_subset$normalized_effect_screen,
  A549_abexinostat_effect_sizes_subset$normalized_effect_sciplex20
)

# Plot the relationship between effect sizes
ggplot(
  A549_abexinostat_effect_sizes_subset,
  aes(x = normalized_effect_screen, y = normalized_effect_sciplex20)
) +
  geom_point(size = 0.5, stroke = 0) +
  geom_smooth(method = "lm", lwd = 0.25) +
  theme(text = element_text(size = 6)) +
  xlab("Screen pseudo-dose DEGs\nbeta coefficients") +
  ylab("Follow-up Pseudo-dose DEG\nbeta coefficients") +
  monocle_theme_opts() +
  ggsave(
    "supplemental_figure32_A549_correlation_of_beta_coefficients_for_union_of_degs_across_Abexinostat_exposed_cells.png",
    height = 2,
    width = 2,
    dpi = 600,
    units = "in"
  )

png(
  "supplemental_figure31_A549_overlap_between_Abexinostat_degs.png",
  width = 3,
  height = 3.25,
  res = 1000,
  units = "in"
)
VennDiagram::draw.pairwise.venn(
  length((
    A549_abexinostat_effect_sizes_subset %>% filter(q_value_screen < 0.05)
  )$id),
  length((
    A549_abexinostat_effect_sizes_subset %>% filter(q_value_sciplex20 < 0.05)
  )$id),
  length(intersect((
    A549_abexinostat_effect_sizes_subset %>% filter(q_value_screen < 0.05)
  )$id,
  (
    A549_abexinostat_effect_sizes_subset %>% filter(q_value_sciplex20 < 0.05)
  )$id
  )),
  category = c("Screen", "sciPlex20"),
  fill = c('red', 'blue'),
  euler.d = TRUE,
  scaled = FALSE,
  cex = 1,
  cat.cex = 1,
  fontfamily = "Helvetica"
)
dev.off()

# For Pracinostat
A549_screen_pracinostat_effect_sizes <-
  A549_screen_HDACi_diff_test_results[["Pracinostat"]] %>%
  dplyr::select(id, gene_short_name, q_value, normalized_effect)

colnames(A549_screen_pracinostat_effect_sizes) <-
  c("id",
    "gene_short_name",
    "q_value_screen",
    "normalized_effect_screen")

A549_sciplex20_pracinostat_effect_sizes <-
  A549_sciplex20_HDACi_diff_test_results[["Pracinostat"]] %>%
  dplyr::select(id, gene_short_name, q_value, normalized_effect)

colnames(A549_sciplex20_pracinostat_effect_sizes) <-
  c("id",
    "gene_short_name",
    "q_value_sciplex20",
    "normalized_effect_sciplex20")

# Ensure data frames are in the same order although left_join should merge safely
identical(
  A549_screen_pracinostat_effect_sizes$id,
  A549_sciplex20_pracinostat_effect_sizes$id
)

# Generate a joint data frame
A549_pracinostat_effect_sizes <-
  left_join(
    A549_screen_pracinostat_effect_sizes,
    A549_sciplex20_pracinostat_effect_sizes[, c("id", "q_value_sciplex20", "normalized_effect_sciplex20")],
    by = "id"
  )

A549_pracinostat_effect_sizes_subset <-
  A549_pracinostat_effect_sizes %>% filter(id %in% A549_Pracinostat_degs)
A549_pracinostat_effect_sizes_subset <-
  A549_pracinostat_effect_sizes_subset %>% filter(normalized_effect_screen != 0 &
                                                    normalized_effect_sciplex20 != 0)

cor.test(
  A549_pracinostat_effect_sizes_subset$normalized_effect_screen,
  A549_pracinostat_effect_sizes_subset$normalized_effect_sciplex20
)

# Plot the relationship between effect sizes
ggplot(
  A549_pracinostat_effect_sizes_subset,
  aes(x = normalized_effect_screen, y = normalized_effect_sciplex20)
) +
  geom_point(size = 0.5, stroke = 0) +
  geom_smooth(method = "lm", lwd = 0.25) +
  theme(text = element_text(size = 6)) +
  xlab("Screen pseudo-dose DEGs\nbeta coefficients") +
  ylab("Follow-up Pseudo-dose DEG\nbeta coefficients") +
  monocle_theme_opts() +
  ggsave(
    "supplemental_figure32_A549_correlation_of_beta_coefficients_for_union_of_degs_across_Pracinostat_exposed_cells.png",
    height = 2,
    width = 2,
    dpi = 600,
    units = "in"
  )

png(
  "supplemental_figure31_A549_Overlap_between_Pracinostat_degs.png",
  width = 3,
  height = 3.25,
  res = 1000,
  units = "in"
)
VennDiagram::draw.pairwise.venn(
  length((
    A549_pracinostat_effect_sizes_subset %>% filter(q_value_screen < 0.05)
  )$id),
  length((
    A549_pracinostat_effect_sizes_subset %>% filter(q_value_sciplex20 < 0.05)
  )$id),
  length(intersect((
    A549_pracinostat_effect_sizes_subset %>% filter(q_value_screen < 0.05)
  )$id,
  (
    A549_pracinostat_effect_sizes_subset %>% filter(q_value_sciplex20 < 0.05)
  )$id
  )),
  category = c("Screen", "sciPlex20"),
  fill = c('red', 'blue'),
  euler.d = TRUE,
  scaled = FALSE,
  cex = 1,
  cat.cex = 1,
  fontfamily = "Helvetica"
)
dev.off()


# For Abexinostat
MCF7_screen_abexinostat_effect_sizes <-
  MCF7_screen_HDACi_diff_test_results[["Abexinostat"]] %>%
  dplyr::select(id, gene_short_name, q_value, normalized_effect)

colnames(MCF7_screen_abexinostat_effect_sizes) <-
  c("id",
    "gene_short_name",
    "q_value_screen",
    "normalized_effect_screen")

MCF7_sciplex20_abexinostat_effect_sizes <-
  MCF7_sciplex20_HDACi_diff_test_results[["Abexinostat"]] %>%
  dplyr::select(id, gene_short_name, q_value, normalized_effect)

colnames(MCF7_sciplex20_abexinostat_effect_sizes) <-
  c("id",
    "gene_short_name",
    "q_value_sciplex20",
    "normalized_effect_sciplex20")

# Ensure data frames are in the same order although left_join should merge safely
identical(
  MCF7_screen_abexinostat_effect_sizes$id,
  MCF7_sciplex20_abexinostat_effect_sizes$id
)

# Generate a joint data frame
MCF7_abexinostat_effect_sizes <-
  left_join(
    MCF7_screen_abexinostat_effect_sizes,
    MCF7_sciplex20_abexinostat_effect_sizes[, c("id", "q_value_sciplex20", "normalized_effect_sciplex20")],
    by = "id"
  )

MCF7_abexinostat_effect_sizes_subset <-
  MCF7_abexinostat_effect_sizes %>% filter(id %in% MCF7_Abexinostat_degs)
MCF7_abexinostat_effect_sizes_subset <-
  MCF7_abexinostat_effect_sizes_subset %>% filter(normalized_effect_screen != 0 &
                                                    normalized_effect_sciplex20 != 0)

cor.test(
  MCF7_abexinostat_effect_sizes_subset$normalized_effect_screen,
  MCF7_abexinostat_effect_sizes_subset$normalized_effect_sciplex20
)

# Plot the relationship between effect sizes
ggplot(MCF7_abexinostat_effect_sizes_subset,
  aes(x = normalized_effect_screen, 
      y = normalized_effect_sciplex20)) +
  geom_point(size = 0.5, stroke = 0) +
  geom_smooth(method = "lm", lwd = 0.25) +
  theme(text = element_text(size = 6)) +
  xlab("Screen pseudo-dose DEGs\nbeta coefficients") +
  ylab("Follow-up Pseudo-dose DEG\nbeta coefficients") +
  monocle_theme_opts() +
  ggsave(
    "supplemental_figure32_MCF7_correlation_of_beta_coefficients_for_union_of_degs_across_Abexinostat_exposed_cells.png",
    height = 2,
    width = 2,
    dpi = 600,
    units = "in"
  )

png(
  "supplemental_figure31_MCF7_overlap_between_Abexinostat_degs.png",
  width = 3,
  height = 3.25,
  res = 1000,
  units = "in"
)
VennDiagram::draw.pairwise.venn(
  length((MCF7_abexinostat_effect_sizes_subset %>% filter(q_value_screen < 0.05))$id),
  length((MCF7_abexinostat_effect_sizes_subset %>% filter(q_value_sciplex20 < 0.05))$id),
  length(intersect((MCF7_abexinostat_effect_sizes_subset %>% filter(q_value_screen < 0.05))$id,
  (MCF7_abexinostat_effect_sizes_subset %>% filter(q_value_sciplex20 < 0.05))$id)),
  category = c("Screen", "sciPlex20"),
  fill = c('red', 'blue'),
  euler.d = TRUE,
  scaled = FALSE,
  cex = 1,
  cat.cex = 1,
  fontfamily = "Helvetica"
)
dev.off()

# For Pracinostat
MCF7_screen_pracinostat_effect_sizes <-
  MCF7_screen_HDACi_diff_test_results[["Pracinostat"]] %>%
  dplyr::select(id, gene_short_name, q_value, normalized_effect)

colnames(MCF7_screen_pracinostat_effect_sizes) <-
  c("id",
    "gene_short_name",
    "q_value_screen",
    "normalized_effect_screen")

MCF7_sciplex20_pracinostat_effect_sizes <-
  MCF7_sciplex20_HDACi_diff_test_results[["Pracinostat"]] %>%
  dplyr::select(id, gene_short_name, q_value, normalized_effect)

colnames(MCF7_sciplex20_pracinostat_effect_sizes) <-
  c("id",
    "gene_short_name",
    "q_value_sciplex20",
    "normalized_effect_sciplex20")

# Ensure data frames are in the same order although left_join should merge safely
identical(
  MCF7_screen_pracinostat_effect_sizes$id,
  MCF7_sciplex20_pracinostat_effect_sizes$id
)

# Generate a joint data frame
MCF7_pracinostat_effect_sizes <-
  left_join(MCF7_screen_pracinostat_effect_sizes,
    MCF7_sciplex20_pracinostat_effect_sizes[, c("id", "q_value_sciplex20", "normalized_effect_sciplex20")],
    by = "id")

MCF7_pracinostat_effect_sizes_subset <-
  MCF7_pracinostat_effect_sizes %>% filter(id %in% MCF7_Pracinostat_degs)
MCF7_pracinostat_effect_sizes_subset <-
  MCF7_pracinostat_effect_sizes_subset %>% 
  filter(normalized_effect_screen != 0,
         normalized_effect_sciplex20 != 0)

cor.test(
  MCF7_pracinostat_effect_sizes_subset$normalized_effect_screen,
  MCF7_pracinostat_effect_sizes_subset$normalized_effect_sciplex20
)

# Plot the relationship between effect sizes
ggplot(MCF7_pracinostat_effect_sizes_subset,
  aes(x = normalized_effect_screen, 
      y = normalized_effect_sciplex20)) +
  geom_point(size = 0.5, stroke = 0) +
  geom_smooth(method = "lm", lwd = 0.25) +
  theme(text = element_text(size = 6)) +
  xlab("Screen pseudo-dose DEGs\nbeta coefficients") +
  ylab("Follow-up Pseudo-dose DEG\nbeta coefficients") +
  monocle_theme_opts() +
  ggsave("supplemental_figure32_MCF7_correlation_of_beta_coefficients_for_union_of_degs_across_Pracinostat_exposed_cells.png",
    height = 2,
    width = 2,
    dpi = 600,
    units = "in")

png("supplemental_figure31_MCF7_Overlap_between_Pracinostat_degs.png",
  width = 3,
  height = 3.25,
  res = 1000,
  units = "in")
VennDiagram::draw.pairwise.venn(
  length((MCF7_pracinostat_effect_sizes_subset %>% filter(q_value_screen < 0.05))$id),
  length((MCF7_pracinostat_effect_sizes_subset %>% filter(q_value_sciplex20 < 0.05))$id),
  length(intersect((MCF7_pracinostat_effect_sizes_subset %>% filter(q_value_screen < 0.05))$id,
  (MCF7_pracinostat_effect_sizes_subset %>% filter(q_value_sciplex20 < 0.05))$id)),
  category = c("Screen", "sciPlex20"),
  fill = c('red', 'blue'),
  euler.d = TRUE,
  scaled = FALSE,
  cex = 1,
  cat.cex = 1,
  fontfamily = "Helvetica"
)
dev.off()
