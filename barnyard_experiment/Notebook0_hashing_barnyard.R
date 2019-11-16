# devtools::install_github(repo = 'cole-trapnell-lab/monocle3', ref = "d4a9a35")
# Specific commit used in performing all analyses

path_to_github = "~/sci-Plex/"

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(stringr)
  library(tidyr)
  library(parallel)
  library(ggridges)
  library(monocle3)
})


DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size = 1000e7)

setwd(paste0(path_to_github,
             "barnyard_experiment",
             sep = ""))

# Set directory for sciPlex github
bin_directory = paste0(path_to_github,
                       "bin",
                       sep = "")




sciPlex_cds <- readRDS("cds_monocle3.RDS")
dim(sciPlex_cds)

sciPlex_cds = sciPlex_cds[, !is.na(pData(sciPlex_cds)$hash_umis)]
colData(sciPlex_cds)$cell = colData(sciPlex_cds)$Cell

cell_meta_data =
  colData(sciPlex_cds) %>%
  as.data.frame() %>%
  pull(cell) %>%
  stringr::str_split_fixed(pattern = "_", 3)

pData(sciPlex_cds)$pcr_well = cell_meta_data[, 1]
pData(sciPlex_cds)$pcr_col =  as.numeric(substr(pData(sciPlex_cds)$pcr_well, 2, 3))
pData(sciPlex_cds)$rt_well = (cell_meta_data[, 3])
pData(sciPlex_cds)$rt_row = substr(pData(sciPlex_cds)$rt_well, 1, 1)
pData(sciPlex_cds)$rt_col = substr(pData(sciPlex_cds)$rt_well, 2, 3)

hash_meta_data = stringr::str_split_fixed(pData(sciPlex_cds)$top_oligo, "_", 2)

pData(sciPlex_cds)$cell_type = hash_meta_data[, 1]
pData(sciPlex_cds)$hash_well = hash_meta_data[, 2]

fresh = c("A", "B", "C", "D")
# frozen = c("E","F","G","H")
pData(sciPlex_cds)$fresh_frozen = ifelse(pData(sciPlex_cds)$rt_row %in% fresh, "fresh", "frozen")

cells.sorted = seq(15, 70, 5)

cells_per_well =
  data.frame(pcr_col = as.numeric(seq(1, 12, 1)),
             cells_sorted = cells.sorted)

colData(sciPlex_cds)$cells_sorted  =
  left_join(colData(sciPlex_cds) %>%
              as.data.frame(),
            cells_per_well, by = "pcr_col") %>%
  pull(cells_sorted)

human.genes = rownames(sciPlex_cds)[grepl("^ENSG", rowData(sciPlex_cds)$id)]
mouse.genes = rownames(sciPlex_cds)[grepl("^ENSMUSG", rowData(sciPlex_cds)$id)]
pData(sciPlex_cds)$n_human_umi = Matrix::colSums(exprs(sciPlex_cds)[rownames(sciPlex_cds) %in% human.genes, ])
pData(sciPlex_cds)$n_mouse_umi = Matrix::colSums(exprs(sciPlex_cds)[rownames(sciPlex_cds) %in% mouse.genes, ])

pData(sciPlex_cds)$h_m_ratio  =
  pData(sciPlex_cds) %>%
  as.data.frame() %>%
  dplyr::mutate(h_m_ratio = n_mouse_umi / (n_mouse_umi + n_human_umi)) %>%
  pull(h_m_ratio)

pData(sciPlex_cds)$h_m_transcriptome =
  sapply(pData(sciPlex_cds)$h_m_ratio, function(x) {
    if (x >= 0.85)
      return("Mouse")
    if (x < 0.15)
      return("Human")
    return("Doublet")
  })

colData(sciPlex_cds) %>%
  as.data.frame() %>%
  ggplot(, aes(x = h_m_ratio)) +
  geom_histogram() +
  geom_vline(xintercept = 0.1) +
  geom_vline(xintercept = 0.9)

colData(sciPlex_cds) %>%
  as.data.frame() %>%
  ggplot()  +
  geom_point(aes(x = n_human_umi, y = n_mouse_umi),
             size = .75,
             stroke = 0) +
  scale_color_manual(values = c("TRUE" = "dimgrey", "FALSE" = "black")) +
  monocle_theme_opts() +
  ggsave("Barnyard_before.png", height = 6, width = 6)

pdf("supplemental_figure1_enrichment_ratio.pdf",
    width = 1.5,
    height = 1.5)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_histogram(aes(x = log10(top_to_second_best_ratio))) +
  monocle_theme_opts() +
  geom_vline(xintercept = log10(15),
             color = "#FF0000",
             size = .5) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("log10(Enrichment Ratio)") +
  ylab("Frequency")
dev.off()

pdf("supplemental_figure1_hashes_per_cell.pdf",
    width = 1.5,
    height = 1.5)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_histogram(aes(x = log10(hash_umis))) +
  monocle_theme_opts() +
  geom_vline(xintercept = log10(10),
             color = "#FF0000",
             size = .5) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("log10(Hash UMIs)") +
  ylab("Frequency")
dev.off()

pdf("figure1_unfiltered_barnyard.pdf",
    width = 1.5,
    height = 1.5)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  ggplot()  +
  geom_point(
    aes(x = n_human_umi, y = n_mouse_umi, color = top_to_second_best_ratio <= 15),
    size = .4,
    stroke = 0
  ) +
  scale_color_manual(
    name = "Detected\nDoublet",
    breaks = c("FALSE", "TRUE"),
    labels = c("Single", "Doublet"),
    values = c("TRUE" = "red", "FALSE" = "black")
  ) +
  monocle3:::monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("Human RNA UMIs") +
  ylab("Mouse RNA UMIs") +
  scale_y_continuous(breaks = c(0, 10000, 20000, 30000),
                     limits = c(0, 25000)) +
  scale_x_continuous(breaks = c(0, 10000, 20000, 30000),
                     limits = c(0, 30000))

dev.off()



pData(sciPlex_cds)$passing_hash_filter =  pData(sciPlex_cds)$top_to_second_best_ratio > 15

cds_filtered = sciPlex_cds[, pData(sciPlex_cds)$passing_hash_filter]

pdf(
  "supplemental_figure1_filtered_barnyard_figure1D.pdf",
  width = 1.5,
  height = 1.5
)
pData(cds_filtered) %>%
  as.data.frame() %>%
  ggplot()  +
  geom_point(aes(x = n_human_umi, y = n_mouse_umi),
             size = .4,
             stroke = 0) +
  monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("Human RNA UMIs") +
  ylab("Mouse RNA UMIs") +
  scale_y_continuous(breaks = c(0, 10000, 20000, 30000),
                     limits = c(0, 25000)) +
  scale_x_continuous(breaks = c(0, 10000, 20000, 30000),
                     limits = c(0, 30000))

dev.off()

pdf(
  "supplemental_figure1_filtered_vs_unfiltered_UMIs.pdf",
  width = 1.5,
  height = 1.5
)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(h_m_transcriptome != "Doublet") %>%
  ggplot() +
  geom_boxplot(
    aes(x = h_m_transcriptome, y = log10(n.umi), fill =  passing_hash_filter),
    size = .25,
    outlier.size = 1,
    outlier.stroke = 0
  ) +
  scale_fill_manual(
    name = "Passed Filter",
    breaks = c("FALSE", "TRUE"),
    labels = c("Single", "Doublet"),
    values = c("TRUE" = "black", "FALSE" = "red")
  ) +
  monocle_theme_opts()  +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("Transcriptome ID") +
  ylab("log10(RNA UMIs)") +
  scale_y_continuous(breaks = c(3, 3.5, 4, 4.5))

dev.off()



# Figure 1E

pdf("figure1_freshfrozen_UMI_Fig1E.pdf",
    width = 1.5,
    height = 1.5)
pData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter, h_m_transcriptome != "Doublet") %>%
  ggplot() +
  geom_boxplot(
    aes(x = h_m_transcriptome, y = log10(n.umi), fill =  fresh_frozen),
    size = .25,
    outlier.size = 1,
    outlier.stroke = 0
  ) +
  scale_fill_manual(
    name = "Passed Filter",
    breaks = c("fresh", "frozen"),
    labels = c("Fresh", "Frozen"),
    values = c("fresh" = "#ef8a62", "frozen" = "#67a9cf")
  ) +
  monocle_theme_opts()  +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("Transcriptome ID") +
  ylab("log10(RNA UMIs)") +
  scale_y_continuous(breaks = c(3, 3.5, 4, 4.5))
dev.off()



pdf(
  "supplemental_figure1_hashes_per_cellType.pdf",
  width = 1.5,
  height = 2
)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter) %>%
  ggplot() +
  geom_boxplot(aes(x = cell_type, y = log10(hash_umis)),
               outlier.size = 1,
               outlier.stroke = 0) +
  monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0,
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  ylab("log10(Hash UMIs)") +
  xlab("")
dev.off()

# This plot shows the hash label is specific for a certain cell type as indicated by the transcriptome
# Supplemental Figure 1D
HEK293_cols = seq(1:6)
NIH3T3_cols = seq(7:12)

pdf("figure1_hash_concordance_Fig1C.pdf",
    width = 1.5,
    height = 1.5)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter) %>%
  mutate(hash_col = as.numeric(substr(hash_well, 2, 3))) %>%
  mutate(human_hash = hash_col %in%  HEK293_cols) %>%
  ggplot() +
  geom_point(
    aes(x = n_human_umi, y = n_mouse_umi, color = human_hash),
    size = .4,
    stroke = 0
  ) +
  monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  scale_color_manual(values = c("TRUE" = "#af8dc3", "FALSE" = "#7fbf7b"),
                     guide = 'none') +
  xlab("Human RNA UMIs") +
  ylab("Mouse RNA UMIs") +
  scale_y_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 30000),
                     limits = c(0, 15000)) +
  scale_x_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 30000),
                     limits = c(0, 20000))
dev.off()



concordant_hash_transcriptome =
  colData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter) %>%
  mutate(hash_col = as.numeric(substr(hash_well, 2, 3))) %>%
  mutate(human_hash = hash_col %in%  HEK293_cols) %>%
  mutate(
    correct_call_hash = (human_hash &
                           h_m_transcriptome == "Human") |
      (!human_hash & h_m_transcriptome == "Mouse")
  ) %>%
  pull(correct_call_hash) %>% sum()

denominator_concordant_hash_transcriptome =
  pData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter) %>%
  pull(cell) %>%
  length()


#concordant_hash_transcriptome_ff =
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  group_by(fresh_frozen) %>%
  filter(passing_hash_filter) %>%
  mutate(hash_col = as.numeric(substr(hash_well, 2, 3))) %>%
  mutate(human_hash = hash_col %in%  HEK293_cols) %>%
  mutate(
    correct_call_hash = (human_hash &
                           h_m_transcriptome == "Human") |
      (!human_hash & h_m_transcriptome == "Mouse")
  ) %>%
  add_tally() %>%
  mutate(rate = sum(correct_call_hash) / n) %>%
  dplyr::select(fresh_frozen, rate) %>%
  distinct()




denominator_concordant_hash_transcriptome =
  colData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter) %>%
  pull(cell) %>%
  length()


print(concordant_hash_transcriptome / denominator_concordant_hash_transcriptome)

colData(cds_filtered) %>%
  as.data.frame() %>%
  filter(h_m_transcriptome == "Doublet")

#################################################################
## Looking at the theoretical doublet rate
#################################################################



doublet_rate_by_cells_sorted <-
  sapply(cells.sorted, function(n.cells.sorted) {
    rt_barcodes = 96
    zero_cells = rt_barcodes * (1 - (1 / rt_barcodes)) ^ n.cells.sorted
    single_cells = n.cells.sorted * (1 - (1 / rt_barcodes)) ^ (n.cells.sorted  - 1)
    
    
    (rt_barcodes - (zero_cells + single_cells)) / n.cells.sorted
    
  })

doublet_rate_by_cells_sorted =
  data.frame(cells.sorted = cells.sorted,
             doublet.rate = doublet_rate_by_cells_sorted)

doublet_rate_by_cells_sorted =
  doublet_rate_by_cells_sorted  %>%
  mutate(
    doublets = cells.sorted * doublet.rate,
    singlets = cells.sorted - doublets,
    singlet.rate = 1 - doublet.rate
  )
rate_df_long =
  doublet_rate_by_cells_sorted %>%
  dplyr::select(everything(),-singlets,-doublets) %>%
  gather(rate.type, rate, doublet.rate, singlet.rate)


pdf(
  "supplemental_figure1_theoretical_proprotion_of_doublets_by_cells_sorted.pdf",
  width = 1.5,
  height = 1.5
)
ggplot(rate_df_long,
       aes(
         x = as.factor(cells.sorted),
         y = rate,
         fill = factor(rate.type, levels = c("singlet.rate", "doublet.rate"))
       )) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    "Cluster Identity",
    labels = c("singlet.rate" = "Singlet", "doublet.rate" = "Multiplet"),
    values = c("singlet.rate" = "black", "doublet.rate" = "red")
  ) +
  theme(
    legend.position = "none",
    legend.text. = ,
    text = element_text(size = 6),
    legend.key.width = unit(0.4, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.text = element_text(size = 6)
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  xlab("Sorted nuclei per well") +
  ylab("Proportion of nuclei") +
  coord_flip() +
  geom_vline(xintercept = 0.25) +
  monocle_theme_opts()
dev.off()

empirical_doublets  =
  colData(sciPlex_cds) %>%
  as.data.frame() %>%
  group_by(cells_sorted, pcr_well) %>%
  dplyr::add_count() %>%
  dplyr::mutate(number_of_doublets = n - sum(passing_hash_filter)) %>%
  dplyr::select(pcr_well, cells_sorted, number_of_doublets) %>%
  dplyr::distinct()


theoretical_doublets =
  rate_df_long %>%
  filter(rate.type == "doublet.rate") %>%
  dplyr::mutate(theoretical_doublets = cells.sorted * rate) %>%
  dplyr::select(cells_sorted = cells.sorted, theoretical_doublets)

cells_per_well_25 = colData(sciPlex_cds) %>% as.data.frame() %>% filter(cells_sorted == 25)
print(paste0(
  "25 cells/well percent passing filter: ",
  100 * sum(cells_per_well_25$passing_hash_filter) / length(cells_per_well_25$passing_hash_filter)
))

cells_per_well_70 = colData(sciPlex_cds) %>% as.data.frame() %>% filter(cells_sorted == 70)
print(paste0("70 cells/well not passing filter: ",
             100 * (
               1 -  sum(cells_per_well_70$passing_hash_filter) / length(cells_per_well_70$passing_hash_filter)
             )))

print(paste0(
  "Percent doublets removed: ",
  100 *
    sum(
      pData(sciPlex_cds)$h_m_transcriptome == "Doublet" &
        !pData(sciPlex_cds)$passing_hash_filter
    ) /
    sum(pData(sciPlex_cds)$h_m_transcriptome == "Doublet")
))

pdf(
  "supplemental_figure1_doublets_by_cells_sorted.pdf",
  width = 1.5,
  height = 1.5
)
left_join(empirical_doublets, theoretical_doublets, by = "cells_sorted") %>%
  ggplot(aes(x = as.factor(cells_sorted))) +
  geom_point(aes(y = number_of_doublets),
             size = .75,
             stroke = 0) + stat_summary(
               aes(y = number_of_doublets),
               fun.y = median,
               fun.ymin = median,
               fun.ymax = median,
               geom = "crossbar",
               width = 0.5,
               fatten = .5,
               colour = "deepskyblue1"
             ) +
  geom_crossbar(
    aes(ymin = theoretical_doublets , y = theoretical_doublets, ymax = theoretical_doublets),
    colour = "red",
    width = 0.5,
    fatten = .5
  ) +
  xlab("Sorted nuclei per well") +
  ylab("Number of Doublets") +
  monocle_theme_opts() +
  theme(
    legend.position = "none",
    legend.text. = ,
    text = element_text(size = 6),
    legend.key.width = unit(0.4, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
dev.off()

rate_df_long %>%
  dplyr::mutate(theoretical_doublets = cells.sorted * rate)

pdf(
  "supplemental_figure1_empirical_proprotion_of_doublets_by_cells_sorted.pdf",
  width = 1.5,
  height = 1.5
)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  group_by(cells_sorted) %>%
  dplyr::add_count() %>%
  dplyr::mutate(proportion_singlet = sum(passing_hash_filter) / n) %>%
  dplyr::select(cells_sorted, proportion_singlet) %>%
  dplyr::distinct() %>%
  dplyr::mutate(proportion_doublet = 1 - proportion_singlet) %>%
  gather(rate_type, rate, proportion_doublet:proportion_singlet) %>%
  ggplot() +
  geom_bar(aes(
    x = as.factor(cells_sorted),
    y = rate,
    fill = factor(
      rate_type,
      levels = c("proportion_singlet", "proportion_doublet")
    )
  ),
  position = "stack",
  stat = "identity") +
  scale_fill_manual(
    "Cluster Identity",
    values = c(
      "proportion_singlet" = "black",
      "proportion_doublet" = "red"
    )
  ) +
  theme(
    legend.position = "none",
    legend.text. = ,
    text = element_text(size = 6),
    legend.key.width = unit(0.4, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.text = element_text(size = 6)
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  xlab("Sorted nuclei per well") +
  ylab("Proportion of nuclei") +
  coord_flip() +
  monocle_theme_opts()
dev.off()


hash_umis_per_cell = read.delim(
  file = "hashUMIs.per.cell",
  header = F,
  col.names = c("sample", "cell", "hash_counts")
)
hash_umis_per_cell =
  hash_umis_per_cell %>%
  dplyr::select(dplyr::everything(),-sample) %>%
  mutate(cell = str_trim(cell, side = c("both")))

rna_umis_per_cell = read.delim(
  file = "UMIs.per.cell.barcode",
  header = F,
  col.names = c("sample", "cell", "rna_counts")
)
rna_umis_per_cell =
  rna_umis_per_cell %>%
  dplyr::select(dplyr::everything(),-sample)

rna_umis_per_cell$cell = str_trim(rna_umis_per_cell$cell, side = c("both"))

rna_hash_union =  left_join(rna_umis_per_cell, hash_umis_per_cell, by = "cell")
rna_hash_union = left_join(rna_hash_union, colData(sciPlex_cds) %>% as.data.frame(), by = "cell")

rna_hash_union$passing_hash_filter = replace_na(rna_hash_union$passing_hash_filter, replace = "low RNA")

pdf("supplemental_figure1_hashes_rna.pdf",
    width = 1.5,
    height = 1.5)
ggplot(rna_hash_union) +
  geom_point(
    aes(
      x = log10(hash_counts) ,
      y = log10(rna_counts),
      color = passing_hash_filter
    ),
    size = .25,
    stroke = 0
  ) +
  monocle_theme_opts() +
  theme(
    legend.position = "none" ,
    text = element_text(size = 6),
    legend.key.width = unit(0.4, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.text = element_text(size = 6)
  ) +
  scale_color_brewer(palette = "Set1") +
  ylab("Log10(RNA UMIs)") +
  xlab("Log10(Hash UMIs)")
dev.off()



pData(sciPlex_cds) %>% filter(passing_hash_filter) %>%
  filter(h_m_transcriptome == "Doublet")

cells_per_hash_well =
  pData(sciPlex_cds) %>%
  as.data.frame() %>%
  mutate(hash_row = substr(hash_well, 1, 1),
         hash_col = substr(hash_well, 2, 3)) %>%
  filter(passing_hash_filter) %>%
  group_by(hash_row, hash_col, cell_type) %>%
  dplyr::summarise(count_per_well = n())



pdf(
  "supplemental_figure1_cells_per_hash_well.pdf",
  width = 3,
  height = 1.5
)
ggplot(cells_per_hash_well) +
  geom_tile(
    aes(
      x = as.factor(hash_row) ,
      y = as.factor(hash_col),
      fill = count_per_well,
      color = cell_type
    ),
    height = .85,
    width = .85,
    size = .35
  ) +
  scale_fill_viridis_c(name = "", limits = c(0, 50)) +
  scale_color_manual(values = c("HEK293T" = "grey80", "NIH3T3" = "black"),
                     guide = 'none') +
  theme(
    legend.position = "right" ,
    text = element_text(size = 6),
    legend.key.width = unit(0.4, "line"),
    legend.key.height = unit(0.6, "line"),
    axis.text = element_text(size = 6)
  ) +
  coord_flip() +
  monocle_theme_opts() +
  ylab("Column") +
  xlab("Row")
dev.off()





pdf("supplemental_figure1_per_hash_well.pdf",
    width = .75,
    height = 1.75)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter) %>%
  group_by(cell_type, hash_well) %>%
  add_tally() %>%
  ggplot() +
  geom_boxplot(
    aes(
      x = as.factor(cell_type) ,
      y = n,
      fill = cell_type
    ),
    size = .25,
    outlier.size = 1,
    outlier.stroke = 0
  ) +
  scale_fill_manual(values = c("HEK293T" = "grey80", "NIH3T3" = "black"),
                    guide = 'none') +
  theme(
    legend.position = "right" ,
    text = element_text(size = 6),
    legend.key.width = unit(0.4, "line"),
    legend.key.height = unit(0.6, "line"),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  monocle_theme_opts()  +
  xlab("") +
  ylab("Cells Recovered")
dev.off()





counts_per_gene =
  pData(cds_filtered) %>%
  as.data.frame() %>%
  group_by(cell_type, fresh_frozen) %>%
  nest() %>%
  mutate(counts_per_gene = purrr::map(
    data,
    .f = function(pdata_subset, cds_filtered) {
      cds_subset = cds_filtered[, pdata_subset$cell]
      totals = Matrix::rowSums(exprs(cds_subset))
      norm_exprs_mat = Matrix::t(Matrix::t(counts(cds_subset)) / size_factors(cds_subset))
      totals = Matrix::rowSums(norm_exprs_mat)
      tibble(gene_id = rownames(exprs(cds_subset)),
             totals = totals)
    },
    cds_filtered
  ))

counts_per_gene =
  counts_per_gene %>%
  unnest(counts_per_gene)

cor_matrix_3T3 =
  counts_per_gene %>%
  filter(cell_type == "NIH3T3") %>%
  spread(key = fresh_frozen, value = totals, fill = 0)

cor(
  x = log10(cor_matrix_3T3$fresh + 1),
  y = log10(cor_matrix_3T3$frozen + 1),
  method = "pearson"
)

cor_matrix_293T =
  counts_per_gene %>%
  filter(cell_type == "HEK293T") %>%
  spread(key = fresh_frozen, value = totals, fill = 0)

cor(
  x = log10(cor_matrix_293T$fresh + 1),
  y = log10(cor_matrix_293T$frozen + 1),
  method = "pearson"
)


counts_per_gene %>%
  spread(key = fresh_frozen, value = totals, fill = 0)  %>%
  ggplot() +
  geom_point(
    aes(x = log10(fresh + 1), y = log10(frozen + 1))
    ,
    stroke = 0,
    size = .25,
    color = "grey80"
  ) +
  geom_abline(intercept = 0,
              slope = 1,
              size = .25) +
  geom_smooth(
    aes(x = log10(fresh + 1), y = log10(frozen + 1)),
    method = "lm",
    color = "red",
    size = .25
  ) +
  stat_cor(
    aes(x = log10(fresh + 1), y = log10(frozen + 1)),
    size = 2,
    label.y.npc = 0.01,
    color = "grey31",
    method = "pearson"
  ) +
  monocle_theme_opts() +
  facet_wrap( ~ cell_type) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    strip.text.x = element_text(size = 6)
  ) +
  xlab("Fresh log10(UMIs + 1)") +
  ylab("Frozen log10(UMIs + 1)") +
  ggsave(
    "supplemental_figure1_fresh_frozen_correlation.png",
    height = 1.5,
    width = 2.25,
    units = "in"
  )

# Figure S1J

pdf("supplemental_figure1_freshfrozen_UMIs.pdf",
    width = 1.5,
    height = 1.5)
colData(sciPlex_cds) %>%
  as.data.frame() %>% ß
filter(passing_hash_filter, h_m_transcriptome != "Doublet") %>%
  ggplot() +
  geom_boxplot(
    aes(
      x = h_m_transcriptome,
      y = log10(hash_umis),
      fill =  fresh_frozen
    ),
    size = .25,
    outlier.size = 1,
    outlier.stroke = 0
  ) +
  scale_fill_manual(
    name = "Passed Filter",
    breaks = c("fresh", "frozen"),
    labels = c("Fresh", "Frozen"),
    values = c("fresh" = "#ef8a62", "frozen" = "#67a9cf")
  ) +
  monocle_theme_opts()  +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("Transcriptome ID") +
  ylab("log10(Hash UMIs)") +
  scale_y_continuous(breaks = c(1, 1.5, 2, 2.5, 3))
dev.off()

#################################################################



sciPlex_cds <- readRDS("sciPlex1_barnyard.RDS")
dim(sciPlex_cds)

sciPlex_cds = sciPlex_cds[, !is.na(pData(sciPlex_cds)$hash_umis)]
colData(sciPlex_cds)$cell = colData(sciPlex_cds)$Cell

cell_meta_data =
  colData(sciPlex_cds) %>%
  as.data.frame() %>%
  pull(cell) %>%
  stringr::str_split_fixed(pattern = "_", 3)

pData(sciPlex_cds)$pcr_well = cell_meta_data[, 1]
pData(sciPlex_cds)$pcr_col =  as.numeric(substr(pData(sciPlex_cds)$pcr_well, 2, 3))
pData(sciPlex_cds)$rt_well = (cell_meta_data[, 3])
pData(sciPlex_cds)$rt_row = substr(pData(sciPlex_cds)$rt_well, 1, 1)
pData(sciPlex_cds)$rt_col = substr(pData(sciPlex_cds)$rt_well, 2, 3)

hash_meta_data = stringr::str_split_fixed(pData(sciPlex_cds)$top_oligo, "_", 2)

pData(sciPlex_cds)$cell_type = hash_meta_data[, 1]
pData(sciPlex_cds)$hash_well = hash_meta_data[, 2]

fresh = c("A", "B", "C", "D")
# frozen = c("E","F","G","H")
pData(sciPlex_cds)$fresh_frozen = ifelse(pData(sciPlex_cds)$rt_row %in% fresh, "fresh", "frozen")

cells.sorted = seq(15, 70, 5)

cells_per_well =
  data.frame(pcr_col = as.numeric(seq(1, 12, 1)),
             cells_sorted = cells.sorted)

colData(sciPlex_cds)$cells_sorted  =
  left_join(colData(sciPlex_cds) %>%
              as.data.frame(),
            cells_per_well, by = "pcr_col") %>%
  pull(cells_sorted)

human.genes = rownames(sciPlex_cds)[grepl("^ENSG", rowData(sciPlex_cds)$id)]
mouse.genes = rownames(sciPlex_cds)[grepl("^ENSMUSG", rowData(sciPlex_cds)$id)]
pData(sciPlex_cds)$n_human_umi = Matrix::colSums(exprs(sciPlex_cds)[rownames(sciPlex_cds) %in% human.genes, ])
pData(sciPlex_cds)$n_mouse_umi = Matrix::colSums(exprs(sciPlex_cds)[rownames(sciPlex_cds) %in% mouse.genes, ])

pData(sciPlex_cds)$h_m_ratio  =
  pData(sciPlex_cds) %>%
  as.data.frame() %>%
  dplyr::mutate(h_m_ratio = n_mouse_umi / (n_mouse_umi + n_human_umi)) %>%
  pull(h_m_ratio)

pData(sciPlex_cds)$h_m_transcriptome =
  sapply(pData(sciPlex_cds)$h_m_ratio, function(x) {
    if (x >= 0.85)
      return("Mouse")
    if (x < 0.15)
      return("Human")
    return("Doublet")
  })

colData(sciPlex_cds) %>%
  as.data.frame() %>%
  ggplot(, aes(x = h_m_ratio)) +
  geom_histogram() +
  geom_vline(xintercept = 0.1) +
  geom_vline(xintercept = 0.9)

colData(sciPlex_cds) %>%
  as.data.frame() %>%
  ggplot()  +
  geom_point(aes(x = n_human_umi, y = n_mouse_umi),
             size = .75,
             stroke = 0) +
  scale_color_manual(values = c("TRUE" = "dimgrey", "FALSE" = "black")) +
  monocle_theme_opts() +
  ggsave("Barnyard_before.png", height = 6, width = 6)

pdf("supplemental_figure1_enrichment_ratio.pdf",
    width = 1.5,
    height = 1.5)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_histogram(aes(x = log10(top_to_second_best_ratio))) +
  monocle_theme_opts() +
  geom_vline(xintercept = log10(15),
             color = "#FF0000",
             size = .5) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("log10(Enrichment Ratio)") +
  ylab("Frequency")
dev.off()

pdf("supplemental_figure1_hashes_per_cell.pdf",
    width = 1.5,
    height = 1.5)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_histogram(aes(x = log10(hash_umis))) +
  monocle_theme_opts() +
  geom_vline(xintercept = log10(10),
             color = "#FF0000",
             size = .5) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("log10(Hash UMIs)") +
  ylab("Frequency")
dev.off()

pdf("figure1_unfiltered_barnyard.pdf",
    width = 1.5,
    height = 1.5)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  ggplot()  +
  geom_point(
    aes(x = n_human_umi, y = n_mouse_umi, color = top_to_second_best_ratio <= 15),
    size = .4,
    stroke = 0
  ) +
  scale_color_manual(
    name = "Detected\nDoublet",
    breaks = c("FALSE", "TRUE"),
    labels = c("Single", "Doublet"),
    values = c("TRUE" = "red", "FALSE" = "black")
  ) +
  monocle3:::monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("Human RNA UMIs") +
  ylab("Mouse RNA UMIs") +
  scale_y_continuous(breaks = c(0, 10000, 20000, 30000),
                     limits = c(0, 25000)) +
  scale_x_continuous(breaks = c(0, 10000, 20000, 30000),
                     limits = c(0, 30000))

dev.off()



pData(sciPlex_cds)$passing_hash_filter =  pData(sciPlex_cds)$top_to_second_best_ratio > 15

cds_filtered = sciPlex_cds[, pData(sciPlex_cds)$passing_hash_filter]

pdf(
  "supplemental_figure1_filtered_barnyard_figure1D.pdf",
  width = 1.5,
  height = 1.5
)
pData(cds_filtered) %>%
  as.data.frame() %>%
  ggplot()  +
  geom_point(aes(x = n_human_umi, y = n_mouse_umi),
             size = .4,
             stroke = 0) +
  monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("Human RNA UMIs") +
  ylab("Mouse RNA UMIs") +
  scale_y_continuous(breaks = c(0, 10000, 20000, 30000),
                     limits = c(0, 25000)) +
  scale_x_continuous(breaks = c(0, 10000, 20000, 30000),
                     limits = c(0, 30000))

dev.off()

pdf(
  "supplemental_figure1_filtered_vs_unfiltered_UMIs.pdf",
  width = 1.5,
  height = 1.5
)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(h_m_transcriptome != "Doublet") %>%
  ggplot() +
  geom_boxplot(
    aes(x = h_m_transcriptome, y = log10(n.umi), fill =  passing_hash_filter),
    size = .25,
    outlier.size = 1,
    outlier.stroke = 0
  ) +
  scale_fill_manual(
    name = "Passed Filter",
    breaks = c("FALSE", "TRUE"),
    labels = c("Single", "Doublet"),
    values = c("TRUE" = "black", "FALSE" = "red")
  ) +
  monocle_theme_opts()  +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("Transcriptome ID") +
  ylab("log10(RNA UMIs)") +
  scale_y_continuous(breaks = c(3, 3.5, 4, 4.5))

dev.off()



# Figure 1E

pdf("figure1_freshfrozen_UMI_Fig1E.pdf",
    width = 1.5,
    height = 1.5)
pData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter, h_m_transcriptome != "Doublet") %>%
  ggplot() +
  geom_boxplot(
    aes(x = h_m_transcriptome, y = log10(n.umi), fill =  fresh_frozen),
    size = .25,
    outlier.size = 1,
    outlier.stroke = 0
  ) +
  scale_fill_manual(
    name = "Passed Filter",
    breaks = c("fresh", "frozen"),
    labels = c("Fresh", "Frozen"),
    values = c("fresh" = "#ef8a62", "frozen" = "#67a9cf")
  ) +
  monocle_theme_opts()  +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("Transcriptome ID") +
  ylab("log10(RNA UMIs)") +
  scale_y_continuous(breaks = c(3, 3.5, 4, 4.5))
dev.off()



pdf(
  "supplemental_figure1_hashes_per_cellType.pdf",
  width = 1.5,
  height = 2
)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter) %>%
  ggplot() +
  geom_boxplot(aes(x = cell_type, y = log10(hash_umis)),
               outlier.size = 1,
               outlier.stroke = 0) +
  monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0,
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  ylab("log10(Hash UMIs)") +
  xlab("")
dev.off()

# This plot shows the hash label is specific for a certain cell type as indicated by the transcriptome
HEK293_cols = seq(1:6)
NIH3T3_cols = seq(7:12)

pdf("figure1_hash_concordance_Fig1C.pdf",
    width = 1.5,
    height = 1.5)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter) %>%
  mutate(hash_col = as.numeric(substr(hash_well, 2, 3))) %>%
  mutate(human_hash = hash_col %in%  HEK293_cols) %>%
  ggplot() +
  geom_point(
    aes(x = n_human_umi, y = n_mouse_umi, color = human_hash),
    size = .4,
    stroke = 0
  ) +
  monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  scale_color_manual(values = c("TRUE" = "#af8dc3", "FALSE" = "#7fbf7b"),
                     guide = 'none') +
  xlab("Human RNA UMIs") +
  ylab("Mouse RNA UMIs") +
  scale_y_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 30000),
                     limits = c(0, 15000)) +
  scale_x_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 30000),
                     limits = c(0, 20000))
dev.off()



concordant_hash_transcriptome =
  colData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter) %>%
  mutate(hash_col = as.numeric(substr(hash_well, 2, 3))) %>%
  mutate(human_hash = hash_col %in%  HEK293_cols) %>%
  mutate(
    correct_call_hash = (human_hash &
                           h_m_transcriptome == "Human") |
      (!human_hash & h_m_transcriptome == "Mouse")
  ) %>%
  pull(correct_call_hash) %>% sum()

denominator_concordant_hash_transcriptome =
  pData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter) %>%
  pull(cell) %>%
  length()


concordant_hash_transcriptome_ff =
  colData(sciPlex_cds) %>%
  as.data.frame() %>%
  group_by(fresh_frozen) %>%
  filter(passing_hash_filter) %>%
  mutate(hash_col = as.numeric(substr(hash_well, 2, 3))) %>%
  mutate(human_hash = hash_col %in%  HEK293_cols) %>%
  mutate(
    correct_call_hash = (human_hash &
                           h_m_transcriptome == "Human") |
      (!human_hash & h_m_transcriptome == "Mouse")
  ) %>%
  add_tally() %>%
  mutate(rate = sum(correct_call_hash) / n) %>%
  dplyr::select(fresh_frozen, rate) %>%
  distinct()




denominator_concordant_hash_transcriptome =
  colData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter) %>%
  pull(cell) %>%
  length()


print(concordant_hash_transcriptome / denominator_concordant_hash_transcriptome)

colData(cds_filtered) %>%
  as.data.frame() %>%
  filter(h_m_transcriptome == "Doublet")

#################################################################
## Looking at the theoretical doublet rate
#################################################################

doublet_rate_by_cells_sorted <-
  sapply(cells.sorted, function(n.cells.sorted) {
    rt_barcodes = 96
    zero_cells = rt_barcodes * (1 - (1 / rt_barcodes)) ^ n.cells.sorted
    single_cells = n.cells.sorted * (1 - (1 / rt_barcodes)) ^ (n.cells.sorted  - 1)
    
    
    (rt_barcodes - (zero_cells + single_cells)) / n.cells.sorted
    
  })

doublet_rate_by_cells_sorted =
  data.frame(cells.sorted = cells.sorted,
             doublet.rate = doublet_rate_by_cells_sorted)

doublet_rate_by_cells_sorted =
  doublet_rate_by_cells_sorted  %>%
  mutate(
    doublets = cells.sorted * doublet.rate,
    singlets = cells.sorted - doublets,
    singlet.rate = 1 - doublet.rate
  )
rate_df_long =
  doublet_rate_by_cells_sorted %>%
  dplyr::select(everything(),-singlets,-doublets) %>%
  gather(rate.type, rate, doublet.rate, singlet.rate)


pdf(
  "supplemental_figure1_theoretical_proprotion_of_doublets_by_cells_sorted.pdf",
  width = 1.5,
  height = 1.5
)
ggplot(rate_df_long,
       aes(
         x = as.factor(cells.sorted),
         y = rate,
         fill = factor(rate.type, levels = c("singlet.rate", "doublet.rate"))
       )) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    "Cluster Identity",
    labels = c("singlet.rate" = "Singlet", "doublet.rate" = "Multiplet"),
    values = c("singlet.rate" = "black", "doublet.rate" = "red")
  ) +
  theme(
    legend.position = "none",
    legend.text. = ,
    text = element_text(size = 6),
    legend.key.width = unit(0.4, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.text = element_text(size = 6)
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  xlab("Sorted nuclei per well") +
  ylab("Proportion of nuclei") +
  coord_flip() +
  geom_vline(xintercept = 0.25) +
  monocle_theme_opts()
dev.off()

empirical_doublets  =
  colData(sciPlex_cds) %>%
  as.data.frame() %>%
  group_by(cells_sorted, pcr_well) %>%
  dplyr::add_count() %>%
  dplyr::mutate(number_of_doublets = n - sum(passing_hash_filter)) %>%
  dplyr::select(pcr_well, cells_sorted, number_of_doublets) %>%
  dplyr::distinct()


theoretical_doublets =
  rate_df_long %>%
  filter(rate.type == "doublet.rate") %>%
  dplyr::mutate(theoretical_doublets = cells.sorted * rate) %>%
  dplyr::select(cells_sorted = cells.sorted, theoretical_doublets)

cells_per_well_25 = colData(sciPlex_cds) %>% as.data.frame() %>% filter(cells_sorted == 25)
print(paste0(
  "25 cells/well percent passing filter: ",
  100 * sum(cells_per_well_25$passing_hash_filter) / length(cells_per_well_25$passing_hash_filter)
))

cells_per_well_70 = colData(sciPlex_cds) %>% as.data.frame() %>% filter(cells_sorted == 70)
print(paste0("70 cells/well not passing filter: ",
             100 * (
               1 -  sum(cells_per_well_70$passing_hash_filter) / length(cells_per_well_70$passing_hash_filter)
             )))

print(paste0(
  "Percent doublets removed: ",
  100 *
    sum(
      pData(sciPlex_cds)$h_m_transcriptome == "Doublet" &
        !pData(sciPlex_cds)$passing_hash_filter
    ) /
    sum(pData(sciPlex_cds)$h_m_transcriptome == "Doublet")
))

pdf(
  "supplemental_figure1_doublets_by_cells_sorted.pdf",
  width = 1.5,
  height = 1.5
)
left_join(empirical_doublets, theoretical_doublets, by = "cells_sorted") %>%
  ggplot(aes(x = as.factor(cells_sorted))) +
  geom_point(aes(y = number_of_doublets),
             size = .75,
             stroke = 0) + stat_summary(
               aes(y = number_of_doublets),
               fun.y = median,
               fun.ymin = median,
               fun.ymax = median,
               geom = "crossbar",
               width = 0.5,
               fatten = .5,
               colour = "deepskyblue1"
             ) +
  geom_crossbar(
    aes(ymin = theoretical_doublets , y = theoretical_doublets, ymax = theoretical_doublets),
    colour = "red",
    width = 0.5,
    fatten = .5
  ) +
  xlab("Sorted nuclei per well") +
  ylab("Number of Doublets") +
  monocle_theme_opts() +
  theme(
    legend.position = "none",
    legend.text. = ,
    text = element_text(size = 6),
    legend.key.width = unit(0.4, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
dev.off()

rate_df_long %>%
  dplyr::mutate(theoretical_doublets = cells.sorted * rate)

pdf(
  "supplemental_figure1_empirical_proprotion_of_doublets_by_cells_sorted.pdf",
  width = 1.5,
  height = 1.5
)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  group_by(cells_sorted) %>%
  dplyr::add_count() %>%
  dplyr::mutate(proportion_singlet = sum(passing_hash_filter) / n) %>%
  dplyr::select(cells_sorted, proportion_singlet) %>%
  dplyr::distinct() %>%
  dplyr::mutate(proportion_doublet = 1 - proportion_singlet) %>%
  gather(rate_type, rate, proportion_doublet:proportion_singlet) %>%
  ggplot() +
  geom_bar(aes(
    x = as.factor(cells_sorted),
    y = rate,
    fill = factor(
      rate_type,
      levels = c("proportion_singlet", "proportion_doublet")
    )
  ),
  position = "stack",
  stat = "identity") +
  scale_fill_manual(
    "Cluster Identity",
    values = c(
      "proportion_singlet" = "black",
      "proportion_doublet" = "red"
    )
  ) +
  theme(
    legend.position = "none",
    legend.text. = ,
    text = element_text(size = 6),
    legend.key.width = unit(0.4, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.text = element_text(size = 6)
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  xlab("Sorted nuclei per well") +
  ylab("Proportion of nuclei") +
  coord_flip() +
  monocle_theme_opts()
dev.off()


hash_umis_per_cell = read.delim(
  file = "hashUMIs.per.cell",
  header = F,
  col.names = c("sample", "cell", "hash_counts")
)
hash_umis_per_cell =
  hash_umis_per_cell %>%
  dplyr::select(dplyr::everything(),-sample) %>%
  mutate(cell = str_trim(cell, side = c("both")))

rna_umis_per_cell = read.delim(
  file = "UMIs.per.cell.barcode",
  header = F,
  col.names = c("sample", "cell", "rna_counts")
)
rna_umis_per_cell =
  rna_umis_per_cell %>%
  dplyr::select(dplyr::everything(),-sample)

rna_umis_per_cell$cell = str_trim(rna_umis_per_cell$cell, side = c("both"))

rna_hash_union =  left_join(rna_umis_per_cell, hash_umis_per_cell, by = "cell")
rna_hash_union = left_join(rna_hash_union, colData(sciPlex_cds) %>% as.data.frame(), by = "cell")

rna_hash_union$passing_hash_filter = replace_na(rna_hash_union$passing_hash_filter, replace = "low RNA")

pdf("supplemental_figure1_hashes_rna.pdf",
    width = 1.5,
    height = 1.5)
ggplot(rna_hash_union) +
  geom_point(
    aes(
      x = log10(hash_counts) ,
      y = log10(rna_counts),
      color = passing_hash_filter
    ),
    size = .25,
    stroke = 0
  ) +
  monocle_theme_opts() +
  theme(
    legend.position = "none" ,
    text = element_text(size = 6),
    legend.key.width = unit(0.4, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.text = element_text(size = 6)
  ) +
  scale_color_brewer(palette = "Set1") +
  ylab("Log10(RNA UMIs)") +
  xlab("Log10(Hash UMIs)")
dev.off()



pData(sciPlex_cds) %>% filter(passing_hash_filter) %>%
  filter(h_m_transcriptome == "Doublet")

cells_per_hash_well =
  pData(sciPlex_cds) %>%
  as.data.frame() %>%
  mutate(hash_row = substr(hash_well, 1, 1),
         hash_col = substr(hash_well, 2, 3)) %>%
  filter(passing_hash_filter) %>%
  group_by(hash_row, hash_col, cell_type) %>%
  dplyr::summarise(count_per_well = n())



pdf(
  "supplemental_figure1_cells_per_hash_well.pdf",
  width = 3,
  height = 1.5
)
ggplot(cells_per_hash_well) +
  geom_tile(
    aes(
      x = as.factor(hash_row) ,
      y = as.factor(hash_col),
      fill = count_per_well,
      color = cell_type
    ),
    height = .85,
    width = .85,
    size = .35
  ) +
  scale_fill_viridis_c(name = "", limits = c(0, 50)) +
  scale_color_manual(values = c("HEK293T" = "grey80", "NIH3T3" = "black"),
                     guide = 'none') +
  theme(
    legend.position = "right" ,
    text = element_text(size = 6),
    legend.key.width = unit(0.4, "line"),
    legend.key.height = unit(0.6, "line"),
    axis.text = element_text(size = 6)
  ) +
  coord_flip() +
  monocle_theme_opts() +
  ylab("Column") +
  xlab("Row")
dev.off()





pdf("supplemental_figure1_per_hash_well.pdf",
    width = .75,
    height = 1.75)
colData(sciPlex_cds) %>%
  as.data.frame() %>%
  filter(passing_hash_filter) %>%
  group_by(cell_type, hash_well) %>%
  add_tally() %>%
  ggplot() +
  geom_boxplot(
    aes(
      x = as.factor(cell_type) ,
      y = n,
      fill = cell_type
    ),
    size = .25,
    outlier.size = 1,
    outlier.stroke = 0
  ) +
  scale_fill_manual(values = c("HEK293T" = "grey80", "NIH3T3" = "black"),
                    guide = 'none') +
  theme(
    legend.position = "right" ,
    text = element_text(size = 6),
    legend.key.width = unit(0.4, "line"),
    legend.key.height = unit(0.6, "line"),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  monocle_theme_opts()  +
  xlab("") +
  ylab("Cells Recovered")
dev.off()


counts_per_gene =
  pData(cds_filtered) %>%
  as.data.frame() %>%
  group_by(cell_type, fresh_frozen) %>%
  nest() %>%
  mutate(counts_per_gene = purrr::map(
    data,
    .f = function(pdata_subset, cds_filtered) {
      cds_subset = cds_filtered[, pdata_subset$cell]
      totals = Matrix::rowSums(exprs(cds_subset))
      norm_exprs_mat = Matrix::t(Matrix::t(counts(cds_subset)) / size_factors(cds_subset))
      totals = Matrix::rowSums(norm_exprs_mat)
      tibble(gene_id = rownames(exprs(cds_subset)),
             totals = totals)
    },
    cds_filtered
  ))

counts_per_gene =
  counts_per_gene %>%
  unnest(counts_per_gene)

cor_matrix_3T3 =
  counts_per_gene %>%
  filter(cell_type == "NIH3T3") %>%
  spread(key = fresh_frozen, value = totals, fill = 0)

cor(
  x = log10(cor_matrix_3T3$fresh + 1),
  y = log10(cor_matrix_3T3$frozen + 1),
  method = "pearson"
)

cor_matrix_293T =
  counts_per_gene %>%
  filter(cell_type == "HEK293T") %>%
  spread(key = fresh_frozen, value = totals, fill = 0)

cor(
  x = log10(cor_matrix_293T$fresh + 1),
  y = log10(cor_matrix_293T$frozen + 1),
  method = "pearson"
)


counts_per_gene %>%
  spread(key = fresh_frozen, value = totals, fill = 0)  %>%
  ggplot() +
  geom_point(
    aes(x = log10(fresh + 1), y = log10(frozen + 1))
    ,
    stroke = 0,
    size = .25,
    color = "grey80"
  ) +
  geom_abline(intercept = 0,
              slope = 1,
              size = .25) +
  geom_smooth(
    aes(x = log10(fresh + 1), y = log10(frozen + 1)),
    method = "lm",
    color = "red",
    size = .25
  ) +
  stat_cor(
    aes(x = log10(fresh + 1), y = log10(frozen + 1)),
    size = 2,
    label.y.npc = 0.01,
    color = "grey31",
    method = "pearson"
  ) +
  monocle_theme_opts() +
  facet_wrap( ~ cell_type) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    strip.text.x = element_text(size = 6)
  ) +
  xlab("Fresh log10(UMIs + 1)") +
  ylab("Frozen log10(UMIs + 1)") +
  ggsave(
    "supplemental_figure1_fresh_frozen_correlation.png",
    height = 1.5,
    width = 2.25,
    units = "in"
  )

# Figure S1J

pdf("supplemental_figure1_freshfrozen_UMIs.pdf",
    width = 1.5,
    height = 1.5)
colData(sciPlex_cds) %>%
  as.data.frame() %>% ß
filter(passing_hash_filter, h_m_transcriptome != "Doublet") %>%
  ggplot() +
  geom_boxplot(
    aes(
      x = h_m_transcriptome,
      y = log10(hash_umis),
      fill =  fresh_frozen
    ),
    size = .25,
    outlier.size = 1,
    outlier.stroke = 0
  ) +
  scale_fill_manual(
    name = "Passed Filter",
    breaks = c("fresh", "frozen"),
    labels = c("Fresh", "Frozen"),
    values = c("fresh" = "#ef8a62", "frozen" = "#67a9cf")
  ) +
  monocle_theme_opts()  +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    legend.text.align = 0
  ) +
  xlab("Transcriptome ID") +
  ylab("log10(Hash UMIs)") +
  scale_y_continuous(breaks = c(1, 1.5, 2, 2.5, 3))
dev.off()

#################################################################
#
#################################################################
