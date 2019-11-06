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
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(reshape2)
  library(pheatmap)
  library(viridis)
  library(Matrix)
  library(devtools)
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

# Figure 3A

# Pull in drug annotations
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

row.names(drug_annotations) = drug_annotations$catalog_number


drug_annotations %>%
  filter(pathway_level_1 != "Vehicle") %>%
  group_by(pathway_level_1) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_bar(
    aes(
      x = reorder(pathway_level_1, n),
      y = n,
      fill = pathway_level_1
    ),
    stat = "identity",
    color = "black",
    size = .25
  ) +
  coord_flip() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.ticks = element_line(colour = 'black', size = .1),
    axis.ticks.length = unit(1, "pt"),
    axis.line.y.left = element_blank(),
    axis.line.x.top = element_blank()
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
  ggsave(
    "Figure3_compounds_by_pathway_level_1.pdf",
    width = 2.5,
    height = 1.5,
    unit = "in",
    dpi = 900
  )










# Load processed sciPlex CDS
sciPlex_cds.list = list()
sciPlex_cds.list[["A549"]] = readRDS("A549_24hrs.RDS")
sciPlex_cds.list[["K562"]] = readRDS("K562_24hrs.RDS")
sciPlex_cds.list[["MCF7"]] = readRDS("MCF7_24hrs.RDS")



product_cluster_mat.list <- list()
cluster.enrichment.df = list()

for (cell_line in names(sciPlex_cds.list)) {
  product_cluster_mat.list[[cell_line]] =
    reshape2::acast(
      colData(sciPlex_cds.list[[cell_line]]) %>%
        as.data.frame() %>%
        dplyr::mutate(dummy = 1) %>%
        dplyr::mutate(product_dose = paste0(product_name, "_", dose)) %>%
        dplyr::select(product_dose, Cluster, dummy),
      product_dose ~ Cluster,
      value.var = "dummy",
      fun.aggregate = sum,
      fill = 0
    )
  
  weighted.mat = product_cluster_mat.list[[cell_line]]
  ntc.counts = product_cluster_mat.list[[cell_line]]["Vehicle_0", ]
  
  cluster.enrichment.df[[cell_line]] = do.call(rbind, lapply(rownames(weighted.mat), function(product_dose) {
    do.call(rbind, lapply(1:ncol(weighted.mat), function(Cluster) {
      test = fisher.test(cbind(c(weighted.mat[product_dose, Cluster], sum(
        weighted.mat[product_dose,-Cluster]
      )),
      c(ntc.counts[Cluster], sum(ntc.counts[-Cluster]))))
      
      data.frame(
        product_dose = product_dose,
        Cluster = Cluster,
        odds.ratio = unname(test$estimate),
        p.value = test$p.value
      )
    }))
  }))
  
  cluster.enrichment.df[[cell_line]]$q.value = p.adjust(cluster.enrichment.df[[cell_line]]$p.value, "BH")
  
  cluster.enrichment.df[[cell_line]]$log2.odds = with(cluster.enrichment.df[[cell_line]],
                                                      ifelse(odds.ratio == 0,-5, round(log2(odds.ratio), 2)))
  
  cluster.enrichment.df[[cell_line]]$product_name <-
    sapply(cluster.enrichment.df[[cell_line]]$product_dose, function(x) {
      stringr::str_split(x, pattern = "_")[[1]][1]
    })
  
  cluster.enrichment.df[[cell_line]]$dose <-
    sapply(cluster.enrichment.df[[cell_line]]$product_dose, function(x) {
      stringr::str_split(x, pattern = "_")[[1]][2]
    })
  
}



significant_product_dose_combninations.list <- list()

for (cell_line in names(sciPlex_cds.list)) {
  significant_product_dose_combninations.list[[cell_line]] <-
    (
      cluster.enrichment.df[[cell_line]] %>%
        filter(q.value < 0.01 &
                 log2.odds > 2.5) %>%
        distinct(product_dose)
    )$product_dose
  
  print(length(significant_product_dose_combninations.list[[cell_line]]))
  
}


min_dose_test.list <- list()
min_dose.list <- list()
min_dose_annotation.list <- list()

for (cell_line in names(sciPlex_cds.list)) {
  min_dose_test.list[[cell_line]] <-
    cluster.enrichment.df[[cell_line]] %>%
    group_by(product_name) %>%
    mutate(min_dose = ifelse(q.value < 0.01 &
                               log2.odds > 2.5, dose, NA)) %>%
    mutate(min_dose = min(min_dose, na.rm = TRUE))
  
  min_dose.list[[cell_line]] <- min_dose_test.list[[cell_line]] %>%
    dplyr::select(product_name, min_dose) %>% distinct()
  
  min_dose_annotation.list[[cell_line]] <-
    sapply(sciPlex_cds.list[[cell_line]]$product_name,
           function(x) {
             min_dose <-
               min_dose.list[[cell_line]][min_dose.list[[cell_line]]$product_name == x, ]$min_dose
             return(min_dose)
             
           })
  
  min_dose_annotation.list[[cell_line]] <-
    as.numeric(min_dose_annotation.list[[cell_line]])
  pData(sciPlex_cds.list[[cell_line]])$min_effective_dose <-
    min_dose_annotation.list[[cell_line]]
  pData(sciPlex_cds.list[[cell_line]])$effective_dose <-
    pData(sciPlex_cds.list[[cell_line]])$dose >= pData(sciPlex_cds.list[[cell_line]])$min_effective_dose
  pData(sciPlex_cds.list[[cell_line]])$min_effective_dose <- NULL
  
}


pathway_level_1_colors =
  c(
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

# Figure 3 UMAP
colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(
    data = colData(sciPlex_cds.list[["A549"]])[!(
      colData(sciPlex_cds.list[["A549"]])$product_dose %in%
        significant_product_dose_combninations.list[["A549"]]
    ),
    c("UMAP_1", "UMAP_2")] %>% as.data.frame(),
    size = 0.1,
    stroke = 0,
    aes(x = UMAP_1, y = UMAP_2),
    color = "gray80",
    alpha = 0.1
  ) +
  geom_point(
    data = colData(sciPlex_cds.list[["A549"]])[colData(sciPlex_cds.list[["A549"]])$product_dose %in%
                                                 significant_product_dose_combninations.list[["A549"]], ] %>% as.data.frame(),
    size = 0.1,
    stroke = 0,
    aes(x = UMAP_1, y = UMAP_2, color = pathway_level_1)
  ) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line")
  ) +
  scale_color_manual("Pathway", values = pathway_level_1_colors) +
  scale_alpha(aes(alpha = dose), range = c(0.01, 1)) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1.5)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "figure3_sciPlex_A549_24hr_UMAP_by_pathway_level_1.png",
    width = 3.25,
    height = 2,
    unit = "in",
    dpi = 600
  )



## Chisq_Cutoff_Residual_model_umi_replicate
colData(sciPlex_cds.list[["K562"]]) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(
    data = as.data.frame(colData(sciPlex_cds.list[["K562"]])[!(
      colData(sciPlex_cds.list[["K562"]])$product_dose %in%
        significant_product_dose_combninations.list[["K562"]]
    ),
    c("UMAP_1", "UMAP_2")]),
    size = 0.1,
    stroke = 0,
    aes(x = UMAP_1, y = UMAP_2),
    color = "gray80",
    alpha = 0.1
  ) +
  geom_point(
    data = as.data.frame(colData(sciPlex_cds.list[["K562"]])[(
      colData(sciPlex_cds.list[["K562"]])$product_dose %in%
        significant_product_dose_combninations.list[["K562"]]
    ), ]),
    size = 0.1,
    stroke = 0,
    aes(x = UMAP_1, y = UMAP_2, color = pathway_level_1)
  ) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line")
  ) +
  scale_color_manual("Pathway", values = pathway_level_1_colors) +
  scale_alpha(aes(alpha = dose), range = c(0.01, 1)) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1.5)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "figure3_sciPlex_K562_24hr_UMAP_by_pathway_level_1.png",
    width = 3.25,
    height = 2,
    unit = "in",
    dpi = 600
  )




## Chisq_Cutoff_Residual_model_umi_replicate
colData(sciPlex_cds.list[["MCF7"]]) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(
    data = as.data.frame(pData(sciPlex_cds.list[["MCF7"]])[!(
      pData(sciPlex_cds.list[["MCF7"]])$product_dose %in%
        significant_product_dose_combninations.list[["MCF7"]]
    ),
    c("UMAP_1", "UMAP_2")]),
    size = 0.1,
    stroke = 0,
    aes(x = UMAP_1, y = UMAP_2),
    color = "gray80",
    alpha = 0.1
  ) +
  geom_point(
    data = as.data.frame(pData(sciPlex_cds.list[["MCF7"]])[pData(sciPlex_cds.list[["MCF7"]])$product_dose %in%
                                                             significant_product_dose_combninations.list[["MCF7"]], ]),
    size = 0.1,
    stroke = 0,
    aes(x = UMAP_1, y = UMAP_2, color = pathway_level_1)
  ) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line")
  ) +
  scale_color_manual("Pathway", values = pathway_level_1_colors) +
  scale_alpha(aes(alpha = dose), range = c(0.01, 1)) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1.5)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "figure3_sciPlex_MCF7_24hr_UMAP_by_pathway_level_1.png",
    width = 3.25,
    height = 2,
    unit = "in",
    dpi = 600
  )



plot_percent_cells_positive(
  sciPlex_cds.list[["A549"]][fData(sciPlex_cds.list[["A549"]])$gene_short_name %in% c("HSP90AA1"),
                             pData(sciPlex_cds.list[["A549"]])$product_name %in% c("Vehicle",
                                                                                   "Trametinib (GSK1120212)")],
  group_cells_by = "dose_character",
  plot_limits = NULL,
  ncol = 3
) +
  scale_fill_manual(
    "Log(Dose [µM])",
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
    axis.text.x = element_text(
      size = 6,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 6)
  ) +
  xlab("Log(Dose [µM])") +
  ylab("Cells (percent)") +
  ggsave(
    "figure3_HSP90AA1_levels_in_trametinib_exposed_A549_cells.png",
    width = 1,
    height = 1,
    unit = "in",
    dpi = 600
  )


plot_percent_cells_positive(
  sciPlex_cds.list[["K562"]][fData(sciPlex_cds.list[["K562"]])$gene_short_name %in% c("HSP90AA1"),
                             pData(sciPlex_cds.list[["K562"]])$product_name %in% c("Vehicle",
                                                                                   "Trametinib (GSK1120212)")],
  group_cells_by = "dose_character",
  plot_limits = NULL,
  ncol = 3
) +
  scale_fill_manual(
    "Log(Dose [µM])",
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
    axis.text.x = element_text(
      size = 6,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 6)
  ) +
  xlab("Log(Dose [µM])") +
  ylab("Cells (percent)") +
  ggsave(
    "figure3_HSP90AA1_levels_in_trametinib_exposed_K562_cells.png",
    width = 1,
    height = 1,
    unit = "in",
    dpi = 600
  )

plot_percent_cells_positive(
  sciPlex_cds.list[["MCF7"]][fData(sciPlex_cds.list[["MCF7"]])$gene_short_name %in% c("HSP90AA1"),
                             pData(sciPlex_cds.list[["MCF7"]])$product_name %in% c("Vehicle",
                                                                                   "Trametinib (GSK1120212)")],
  group_cells_by = "dose_character",
  plot_limits = NULL,
  ncol = 3
) +
  scale_fill_manual(
    "Log(Dose [µM])",
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
    axis.text.x = element_text(
      size = 6,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 6)
  ) +
  xlab("Log(Dose [µM])") +
  ylab("Cells (percent)") +
  ggsave(
    "figure3_HSP90AA1_levels_in_trametinib_exposed_MCF7_cells.png",
    width = 1,
    height = 1,
    unit = "in",
    dpi = 600
  )


ggplot(as.data.frame((colData(sciPlex_cds.list[["A549"]])))[pData(sciPlex_cds.list[["A549"]])$product_name %in%
                                                              c(
                                                                "Vehicle",
                                                                "Trametinib (GSK1120212)",
                                                                "Alvespimycin (17-DMAG) HCl",
                                                                "Tanespimycin (17-AAG)",
                                                                "Luminespib (AUY-922, NVP-AUY922)"
                                                              ), ],
       aes(
         x = UMAP_1,
         y = UMAP_2,
         color = factor(
           product_name,
           levels = c(
             "Vehicle",
             "Trametinib (GSK1120212)",
             "Alvespimycin (17-DMAG) HCl",
             "Tanespimycin (17-AAG)",
             "Luminespib (AUY-922, NVP-AUY922)"
           )
         )
       )) +
  geom_point(size = 0.15, stroke = 0) +
  scale_color_manual(
    "Treatment",
    values = c(
      "Vehicle" = "grey80",
      "Trametinib (GSK1120212)" = "firebrick1",
      "Alvespimycin (17-DMAG) HCl" = "darkorchid2",
      "Tanespimycin (17-AAG)" = "darkorchid2",
      "Luminespib (AUY-922, NVP-AUY922)" = "darkorchid2"
    )
  ) +
  #geom_density2d() +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 2)))) +
  theme_void() +
  theme(legend.position = "none") +
  ggsave(
    "figure3_sciPlex_A549_24hr_UMAP_by_Trametinib_Residual_model_umi_replicate.png",
    width = 1,
    height = 1,
    unit = "in",
    dpi = 600
  )

ggplot(as.data.frame(colData(sciPlex_cds.list[["K562"]]))[(
  pData(sciPlex_cds.list[["K562"]])$product_name %in%
    c(
      "Vehicle",
      "Trametinib (GSK1120212)",
      "Alvespimycin (17-DMAG) HCl",
      "Tanespimycin (17-AAG)",
      "Luminespib (AUY-922, NVP-AUY922)"
    )
) , ],
aes(
  x = UMAP_1,
  y = UMAP_2,
  color = factor(
    product_name,
    levels = c(
      "Vehicle",
      "Trametinib (GSK1120212)",
      "Alvespimycin (17-DMAG) HCl",
      "Tanespimycin (17-AAG)",
      "Luminespib (AUY-922, NVP-AUY922)"
    )
  )
)) +
  geom_point(size = 0.15, stroke = 0) +
  scale_color_manual(
    "Treatment",
    values = c(
      "Vehicle" = "grey80",
      "Trametinib (GSK1120212)" = "firebrick1",
      "Alvespimycin (17-DMAG) HCl" = "darkorchid2",
      "Tanespimycin (17-AAG)" = "darkorchid2",
      "Luminespib (AUY-922, NVP-AUY922)" = "darkorchid2"
    )
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  theme_void() +
  theme(legend.position = "none") +
  ggsave(
    "figure3_sciPlex_K562_24hr_UMAP_by_Trametinib_Residual_model_umi_replicate.png",
    width = 1,
    height = 1,
    unit = "in",
    dpi = 600
  )

ggplot(as.data.frame(pData(sciPlex_cds.list[["MCF7"]]))[(
  pData(sciPlex_cds.list[["MCF7"]])$product_name %in%
    c(
      "Vehicle",
      "Trametinib (GSK1120212)",
      "Alvespimycin (17-DMAG) HCl",
      "Tanespimycin (17-AAG)",
      "Luminespib (AUY-922, NVP-AUY922)"
    )
) , ],
aes(
  x = UMAP_1,
  y = UMAP_2,
  color = factor(
    product_name,
    levels = c(
      "Vehicle",
      "Trametinib (GSK1120212)",
      "Alvespimycin (17-DMAG) HCl",
      "Tanespimycin (17-AAG)",
      "Luminespib (AUY-922, NVP-AUY922)"
    )
  )
)) +
  geom_point(size = 0.15, stroke = 0) +
  scale_color_manual(
    "Treatment",
    values = c(
      "Vehicle" = "grey80",
      "Trametinib (GSK1120212)" = "firebrick1",
      "Alvespimycin (17-DMAG) HCl" = "darkorchid2",
      "Tanespimycin (17-AAG)" = "darkorchid2",
      "Luminespib (AUY-922, NVP-AUY922)" = "darkorchid2"
    )
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  theme_void() +
  theme(legend.position = "none") +
  ggsave(
    "figure3_sciPlex_MCF7_24hr_UMAP_by_Trametinib_Residual_model_umi_replicate.png",
    width = 1,
    height = 1,
    unit = "in",
    dpi = 600
  )


sciPlex_counts.list <- list()
sciPlex_models.list <- list()
sciPlex_model_stats.list <- list()

for (cell_line in names(sciPlex_cds.list)) {
  # Count the number of cells from each well and normalize per replicate and culture plate
  sciPlex_counts.list[[cell_line]] =
    count_cells(
      sciPlex_cds.list[[cell_line]],
      normalization_terms = c("cell_type", "culture_plate", "replicate")
    )
  
  sciPlex_models.list[[cell_line]] =
    fit_dose_response_models(sciPlex_counts.list[[cell_line]])
  
  sciPlex_model_stats.list[[cell_line]] =
    extract_dose_model_stats(sciPlex_counts.list[[cell_line]], sciPlex_models.list[[cell_line]])
  
}


A549_viability =
  colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  filter(!vehicle) %>%
  dplyr::select(product_name, dose, viability) %>%
  distinct() %>%
  arrange(product_name, dose) %>%
  spread(key = product_name, value = viability)

row.names(A549_viability) <-
  c("A549_10", "A549_100", "A549_1000", "A549_10000")

MCF7_viability =
  colData(sciPlex_cds.list[["MCF7"]]) %>%
  as.data.frame() %>%
  filter(!vehicle) %>%
  dplyr::select(product_name, dose, viability) %>%
  distinct() %>%
  arrange(product_name, dose) %>%
  spread(key = product_name, value = viability)

row.names(MCF7_viability) <-
  c("MCF7_10", "MCF7_100", "MCF7_1000", "MCF7_10000")

K562_viability <-
  colData(sciPlex_cds.list[["K562"]]) %>%
  as.data.frame() %>%
  filter(!vehicle) %>%
  dplyr::select(product_name, dose, viability) %>%
  distinct() %>%
  arrange(product_name, dose) %>%
  spread(key = product_name, value = viability)

row.names(K562_viability) <-
  c("K562_10", "K562_100", "K562_1000", "K562_10000")


A549_dose_dependent_drugs <-
  colnames(A549_viability[, apply(A549_viability < 0.5, 2, any)])
MCF7_dose_dependent_drugs <-
  colnames(MCF7_viability[, apply(MCF7_viability < 0.5, 2, any)])
K562_dose_dependent_drugs <-
  colnames(K562_viability[, apply(K562_viability < 0.5, 2, any)])


dose_dependent_drugs <- Reduce(
  union,
  list(
    A549_dose_dependent_drugs,
    MCF7_dose_dependent_drugs,
    K562_dose_dependent_drugs
  )
)


A549_viability_subset <- A549_viability[, dose_dependent_drugs]
MCF7_viability_subset <- MCF7_viability[, dose_dependent_drugs]
K562_viability_subset <- K562_viability[, dose_dependent_drugs]

dim(A549_viability_subset)
dim(MCF7_viability_subset)
dim(K562_viability_subset)

identical(colnames(A549_viability_subset),
          colnames(MCF7_viability_subset))
identical(colnames(A549_viability_subset),
          colnames(K562_viability_subset))

viability_matrix <- rbind(A549_viability_subset,
                          MCF7_viability_subset,
                          K562_viability_subset)

length(dose_dependent_drugs)

dose_dependent_drug_metadata =
  colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  filter(product_name %in% dose_dependent_drugs) %>%
  dplyr::select(product_name, target, pathway_level_1, pathway_level_2) %>%
  distinct()

dose_dependent_drug_metadata$target <-
  as.character(dose_dependent_drug_metadata$target)
dose_dependent_drug_metadata$pathway_level_1 <-
  as.character(dose_dependent_drug_metadata$pathway_level_1)
dose_dependent_drug_metadata$pathway_level_2 <-
  as.character(dose_dependent_drug_metadata$pathway_level_2)

dose_dependent_drug_metadata$target[dose_dependent_drug_metadata$product_name == "Bosutinib (SKI-606)"] <-
  "Bcr-Abl,Src"

dose_dependent_drug_metadata$pathway_level_2[dose_dependent_drug_metadata$pathway_level_2 == "Tyrosine kinase activity"] <-
  "RTK activity"


dose_dependent_drug_metadata$target[dose_dependent_drug_metadata$product_name %in%
                                      c("Navitoclax (ABT-263)")] <-
  "Bcl2"


dose_dependent_drug_metadata$pathway_level_2[dose_dependent_drug_metadata$product_name %in%
                                               c("Bosutinib (SKI-606)")] <-
  "Abl/Src activity"


sorted_drugs = as.character((dose_dependent_drug_metadata %>% arrange(pathway_level_2))$product_name)

head(sorted_drugs)

ann_colors <-
  list(
    "Pathway" = c(
      "JAK kinase activity" = "aquamarine3",
      "RTK activity" = "lemonchiffon4",
      "Histone deacetylation" = "deepskyblue2",
      "Mitochondria-mediated apoptosis" = "slategray4",
      "Abl/Src activity" = "firebrick1",
      "HSP90 activity" = "brown4",
      "MAPK activity" = "darkgreen",
      "DNA methylation" = "orangered3",
      
      "Cell cycle regulation" = "gold4",
      "Aurora kinase activity" = "darkolivegreen",
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
  )

length(unique(dose_dependent_drug_metadata$pathway_level_2))

ph <- pheatmap(
  viability_matrix[, sorted_drugs],
  clustering_method = "ward.D2",
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = magma(35, direction = -1),
  annotation_col = data.frame(
    row.names = dose_dependent_drug_metadata$product_name,
    Pathway = dose_dependent_drug_metadata$pathway_level_2
  ),
  annotation_colors = ann_colors ,
  treeheight_row = 10,
  treeheight_col = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  file = "figure3_Cell_type_viability_drug.png",
  width = 10,
  height = 10
)

ph <- pheatmap(
  viability_matrix[, sorted_drugs],
  clustering_method = "ward.D2",
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = magma(35, direction = -1),
  annotation_col = data.frame(
    row.names = dose_dependent_drug_metadata$product_name,
    Pathway = dose_dependent_drug_metadata$pathway_level_2
  ),
  annotation_colors = ann_colors ,
  treeheight_row = 10,
  treeheight_col = 10,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  file = "figure3_Cell_type_viability_drug_clustered.png",
  width = 10,
  height = 5
)

pheatmap(
  A549_viability_subset[, sorted_drugs],
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = magma(35, direction = -1),
  annotation_col = data.frame(
    row.names = dose_dependent_drug_metadata$product_name,
    Pathway = dose_dependent_drug_metadata$pathway_level_2
  ),
  annotation_colors = ann_colors ,
  treeheight_row = 10,
  treeheight_col = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  file = "figure3_A549_viability_by_drug.png",
  cutree_rows = 3.25,
  width = 12,
  height = 3.75
)

pheatmap(
  MCF7_viability_subset[, sorted_drugs],
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = magma(35, direction = -1),
  annotation_col = data.frame(
    row.names = dose_dependent_drug_metadata$product_name,
    Pathway = dose_dependent_drug_metadata$pathway_level_2
  ),
  annotation_colors = ann_colors ,
  treeheight_row = 10,
  treeheight_col = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  file = "figure3_MCF7_viability_by_drug.png",
  cutree_rows = 3.25,
  width = 12,
  height = 3.75
)

pheatmap(
  K562_viability_subset[, sorted_drugs],
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = magma(35, direction = -1),
  annotation_col = data.frame(
    row.names = dose_dependent_drug_metadata$product_name,
    Pathway = dose_dependent_drug_metadata$pathway_level_2
  ),
  annotation_colors = ann_colors ,
  treeheight_row = 10,
  treeheight_col = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  file = "figure3_K562_viability_by_drug.png",
  cutree_rows = 3.25,
  width = 12,
  height = 3.75
)


######################################

# Supplemental Figure Cluster Proportions Cell Cycle Regulators

A549_cell_proportion_by_cell_cycle_drugs <-
  colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  filter(pathway_level_1 %in% c("Vehicle", "Cell cycle regulation")) %>%
  dplyr::mutate(pathway = factor(pathway_level_1, levels = c("Vehicle", "Cell cycle regulation"))) %>%
  group_by(pathway, Cluster) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  complete(pathway, Cluster, fill = list(n = 0)) %>%
  group_by(pathway) %>%
  dplyr::mutate(proportion = n / sum(n))

K562_cell_proportion_by_cell_cycle_drugs <-
  colData(sciPlex_cds.list[["K562"]]) %>%
  as.data.frame() %>%
  filter(pathway_level_1 %in% c("Vehicle", "Cell cycle regulation")) %>%
  dplyr::mutate(pathway = factor(pathway_level_1, levels = c("Vehicle", "Cell cycle regulation"))) %>%
  group_by(pathway, Cluster) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  complete(pathway, Cluster, fill = list(n = 0)) %>%
  group_by(pathway) %>%
  dplyr::mutate(proportion = n / sum(n))

MCF7_cell_proportion_by_cell_cycle_drugs <-
  colData(sciPlex_cds.list[["MCF7"]]) %>%
  as.data.frame() %>%
  filter(pathway_level_1 %in% c("Vehicle", "Cell cycle regulation")) %>%
  dplyr::mutate(pathway = factor(pathway_level_1, levels = c("Vehicle", "Cell cycle regulation"))) %>%
  group_by(pathway, Cluster) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  complete(pathway, Cluster, fill = list(n = 0)) %>%
  group_by(pathway) %>%
  dplyr::mutate(proportion = n / sum(n))


A549_clusters_ordered <-
  (
    A549_cell_proportion_by_cell_cycle_drugs %>%
      filter(pathway == "Vehicle") %>%
      arrange(desc(proportion))
  )$Cluster

A549_cell_proportion_by_cell_cycle_drugs$Cluster <-
  factor(A549_cell_proportion_by_cell_cycle_drugs$Cluster,
         levels = A549_clusters_ordered)
ggplot(A549_cell_proportion_by_cell_cycle_drugs,
       aes(x = Cluster, y = proportion)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ factor(pathway,
                       levels = c("Vehicle", "Cell cycle regulation")),
              scales = "free") +
  coord_flip() +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    axis.text.x = element_text(
      size = 4,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 4)
  ) +
  ylab("Proportion of cells") +
  monocle:::monocle_theme_opts() +
  ggsave(
    "A549_distribution_of_cell_cycle_drugs_across_PCA_clusters.png",
    width = 2,
    height = 2
  )

MCF7_clusters_ordered <-
  (
    MCF7_cell_proportion_by_cell_cycle_drugs %>%
      filter(pathway == "Vehicle") %>%
      arrange(desc(proportion))
  )$Cluster

MCF7_cell_proportion_by_cell_cycle_drugs$Cluster <-
  factor(MCF7_cell_proportion_by_cell_cycle_drugs$Cluster,
         levels = MCF7_clusters_ordered)

ggplot(MCF7_cell_proportion_by_cell_cycle_drugs,
       aes(x = Cluster, y = proportion)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ factor(pathway,
                       levels = c("Vehicle", "Cell cycle regulation")),
              scales = "free") +
  coord_flip() +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    axis.text.x = element_text(
      size = 4,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 4)
  ) +
  ylab("Proportion of cells") +
  monocle:::monocle_theme_opts() +
  ggsave(
    "MCF7_distribution_of_cell_cycle_drugs_across_PCA_clusters.png",
    width = 2,
    height = 2
  )

K562_clusters_ordered <-
  (
    K562_cell_proportion_by_cell_cycle_drugs %>%
      filter(pathway == "Vehicle") %>%
      arrange(desc(proportion))
  )$Cluster

K562_cell_proportion_by_cell_cycle_drugs$Cluster <-
  factor(K562_cell_proportion_by_cell_cycle_drugs$Cluster,
         levels = K562_clusters_ordered)

ggplot(K562_cell_proportion_by_cell_cycle_drugs,
       aes(x = Cluster, y = proportion)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ factor(pathway,
                       levels = c("Vehicle", "Cell cycle regulation")),
              scales = "free") +
  coord_flip() +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    axis.text.x = element_text(
      size = 4,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 4)
  ) +
  ylab("Proportion of cells") +
  monocle:::monocle_theme_opts() +
  ggsave(
    "K562_distribution_of_cell_cycle_drugs_across_PCA_clusters.png",
    width = 2,
    height = 2
  )



# Supplemental Figure Cluster Proportions

# Residual_model_umi_replicate
ggplot(
  colData(sciPlex_cds.list[["A549"]])  %>% as.data.frame(),
  aes(x = UMAP_1, y = UMAP_2, color = Cluster)
) +
  geom_point(size = 0.05, stroke = 0) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line")
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure7_sciPlex_A549_24hr_UMAP_by_Louvain_Cluster.png",
    width = 3,
    height = 2.25
  )

# Residual_model_umi_replicate
ggplot(
  colData(sciPlex_cds.list[["K562"]]) %>% as.data.frame(),
  aes(x = UMAP_1, y = UMAP_2, color = Cluster)
) +
  geom_point(size = 0.05, stroke = 0) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line")
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle:::monocle_theme_opts() +
  ggsave(
    "supplemental_figure7_sciPlex_K562_24hr_UMAP_by_Louvain_Cluster.png",
    width = 3,
    height = 2.25
  )

# Residual_model_umi_replicate
ggplot(
  colData(sciPlex_cds.list[["MCF7"]]) %>% as.data.frame(),
  aes(x = UMAP_1, y = UMAP_2, color = Cluster)
) +
  geom_point(size = 0.05, stroke = 0) +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line")
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle_theme_opts() +
  ggsave(
    "supplemental_figure7_sciPlex_MCF7_24hr_UMAP_by_Louvain_Cluster.png",
    width = 3,
    height = 2.25
  )

# Triamcinalone treatment of A549 cells
ggplot(colData(sciPlex_cds.list[["A549"]]) %>% as.data.frame(),
       aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    alpha = 0.5,
    color = "grey80",
    size = 0.05,
    stroke = 0
  ) +
  geom_point(
    data =
      colData(sciPlex_cds.list[["A549"]]) %>%
      as.data.frame() %>%
      filter(catalog_number == "S1628"),
    color = "chartreuse3",
    size = 0.05,
    stroke = 0
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line")
  ) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  monocle_theme_opts() +
  ggsave(
    "supplemental_figure7_sciPlex_A549_24hr_UMAP_by_triamcinalone.png",
    width = 1.5,
    height = 1.5
  )

#############

# Epothilone A by dose in A549 cells
ggplot(colData(sciPlex_cds.list[["A549"]]) %>% as.data.frame(),
       aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "grey80",
             size = 0.05,
             stroke = 0) +
  geom_point(
    data =
      colData(sciPlex_cds.list[["A549"]]) %>%
      as.data.frame() %>%
      filter(catalog_number == "S1297"),
    aes(color = dose_character),
    size = 0.15,
    stroke = 0
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text =  element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_blank()
  ) +
  scale_color_manual(
    "Dose [nM]",
    values = c(
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  ggsave(
    "supplemental_figure7_sciPlex_A549_24hr_UMAP_by_epoA_by_dose.png",
    width = 1.5,
    height = 1.5,
    dpi = 900
  )



# Epothilone A by dose in A549 cells
ggplot(
  colData(sciPlex_cds.list[["A549"]]) %>%
    as.data.frame() %>%
    filter(UMAP_1 > 2,
           UMAP_1 < 4,
           UMAP_2 < -0.5,
           UMAP_2 > -2.5),
  aes(x = UMAP_1, y = UMAP_2)
) +
  geom_point(color = "grey80",
             size = 0.25,
             stroke = 0) +
  geom_point(
    data =
      colData(sciPlex_cds.list[["A549"]]) %>%
      as.data.frame() %>%
      filter(
        catalog_number == "S1297",
        UMAP_1 > 2,
        UMAP_1 < 4,
        UMAP_2 < -0.5,
        UMAP_2 > -2.5
      ),
    aes(color = dose_character),
    size = 0.5,
    stroke = 0
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text =  element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_blank()
  ) +
  scale_color_manual(
    "Dose [nM]",
    values = c(
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  ggsave(
    "supplemental_figure7_sciPlex_A549_24hr_UMAP_by_epoA_sub1.png",
    width = .75,
    height = .75,
    dpi = 900
  )


# Epothilone A by dose in A549 cells
ggplot(
  colData(sciPlex_cds.list[["A549"]]) %>%
    as.data.frame() %>%
    filter(UMAP_1 > 3,
           UMAP_1 < 4,
           UMAP_2 < 2,
           UMAP_2 > 1),
  aes(x = UMAP_1, y = UMAP_2)
) +
  geom_point(color = "grey80",
             size = 0.25,
             stroke = 0) +
  geom_point(
    data =
      colData(sciPlex_cds.list[["A549"]]) %>%
      as.data.frame() %>%
      filter(
        catalog_number == "S1297",
        UMAP_1 > 3,
        UMAP_1 < 4,
        UMAP_2 < 2,
        UMAP_2 > 1
      ),
    aes(color = dose_character),
    size = 0.5,
    stroke = 0
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text =  element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_blank()
  ) +
  scale_color_manual(
    "Dose [nM]",
    values = c(
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  ggsave(
    "supplemental_figure7_sciPlex_A549_24hr_UMAP_by_epoA_sub2.png",
    width = .75 / 2,
    height = .75 / 2,
    dpi = 900
  )


# Epothilone A by dose in A549 cells
ggplot(
  colData(sciPlex_cds.list[["A549"]]) %>%
    as.data.frame() %>%
    filter(UMAP_1 > -6.5,
           UMAP_1 < -4.5,
           UMAP_2 < 1,
           UMAP_2 > -1),
  aes(x = UMAP_1, y = UMAP_2)
) +
  geom_point(color = "grey80",
             size = 0.25,
             stroke = 0) +
  geom_point(
    data =
      colData(sciPlex_cds.list[["A549"]]) %>%
      as.data.frame() %>%
      filter(
        catalog_number == "S1297",
        UMAP_1 > -6.5,
        UMAP_1 < -4.5,
        UMAP_2 < 1,
        UMAP_2 > -1
      ),
    aes(color = dose_character),
    size = 0.5,
    stroke = 0
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text =  element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_blank()
  ) +
  scale_color_manual(
    "Dose [nM]",
    values = c(
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  ggsave(
    "supplemental_figure7_sciPlex_A549_24hr_UMAP_by_epoA_sub3.png",
    width = .75,
    height = .75,
    dpi = 900
  )


# Epothilone B by dose in A549 cells
ggplot(colData(sciPlex_cds.list[["A549"]]) %>% as.data.frame(),
       aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "grey80",
             size = 0.05,
             stroke = 0) +
  geom_point(
    data =
      colData(sciPlex_cds.list[["A549"]]) %>%
      as.data.frame() %>%
      filter(catalog_number == "S1364"),
    aes(color = dose_character),
    size = 0.15,
    stroke = 0
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text =  element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_blank()
  ) +
  scale_color_manual(
    "Dose [nM]",
    values = c(
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  ggsave(
    "supplemental_figure7_sciPlex_A549_24hr_UMAP_by_epoB_by_dose.png",
    width = 1.5,
    height = 1.5,
    dpi = 900
  )

# Epothilone B by dose in A549 cells
ggplot(colData(sciPlex_cds.list[["A549"]]) %>% as.data.frame(),
       aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "grey80",
             size = 0.05,
             stroke = 0) +
  geom_point(
    data =
      colData(sciPlex_cds.list[["A549"]]) %>%
      as.data.frame() %>%
      filter(catalog_number == "S1364"),
    aes(color = dose_character),
    size = 0.15,
    stroke = 0
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text =  element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_blank()
  ) +
  scale_color_manual(
    "Dose [nM]",
    values = c(
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  ggsave(
    "supplemental_figure7_sciPlex_A549_24hr_UMAP_by_epoB_by_dose_legend.png",
    width = 2,
    height = 1.5,
    dpi = 900
  )



# Epothilone A by dose in A549 cells
ggplot(
  colData(sciPlex_cds.list[["A549"]]) %>%
    as.data.frame() %>%
    filter(UMAP_1 > 2,
           UMAP_1 < 4,
           UMAP_2 < -0.5,
           UMAP_2 > -2.5),
  aes(x = UMAP_1, y = UMAP_2)
) +
  geom_point(color = "grey80",
             size = 0.25,
             stroke = 0) +
  geom_point(
    data =
      colData(sciPlex_cds.list[["A549"]]) %>%
      as.data.frame() %>%
      filter(
        catalog_number == "S1364",
        UMAP_1 > 2,
        UMAP_1 < 4,
        UMAP_2 < -0.5,
        UMAP_2 > -2.5
      ),
    aes(color = dose_character),
    size = 0.5,
    stroke = 0
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text =  element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_blank()
  ) +
  scale_color_manual(
    "Dose [nM]",
    values = c(
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  ggsave(
    "supplemental_figure7_sciPlex_A549_24hr_UMAP_by_epoB_sub1.png",
    width = .75,
    height = .75,
    dpi = 900
  )


# Epothilone A by dose in A549 cells
ggplot(
  colData(sciPlex_cds.list[["A549"]]) %>%
    as.data.frame() %>%
    filter(UMAP_1 > 3,
           UMAP_1 < 4,
           UMAP_2 < 2,
           UMAP_2 > 1),
  aes(x = UMAP_1, y = UMAP_2)
) +
  geom_point(color = "grey80",
             size = 0.25,
             stroke = 0) +
  geom_point(
    data =
      colData(sciPlex_cds.list[["A549"]]) %>%
      as.data.frame() %>%
      filter(
        catalog_number == "S1364",
        UMAP_1 > 3,
        UMAP_1 < 4,
        UMAP_2 < 2,
        UMAP_2 > 1
      ),
    aes(color = dose_character),
    size = 0.5,
    stroke = 0
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text =  element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_blank()
  ) +
  scale_color_manual(
    "Dose [nM]",
    values = c(
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  ggsave(
    "supplemental_figure7_sciPlex_A549_24hr_UMAP_by_epoB_sub2.png",
    width = .75 / 2,
    height = .75 / 2,
    dpi = 900
  )


# Epothilone A by dose in A549 cells
ggplot(
  colData(sciPlex_cds.list[["A549"]]) %>%
    as.data.frame() %>%
    filter(UMAP_1 > -6.5,
           UMAP_1 < -4.5,
           UMAP_2 < 1,
           UMAP_2 > -1),
  aes(x = UMAP_1, y = UMAP_2)
) +
  geom_point(color = "grey80",
             size = 0.25,
             stroke = 0) +
  geom_point(
    data =
      colData(sciPlex_cds.list[["A549"]]) %>%
      as.data.frame() %>%
      filter(
        catalog_number == "S1364",
        UMAP_1 > -6.5,
        UMAP_1 < -4.5,
        UMAP_2 < 1,
        UMAP_2 > -1
      ),
    aes(color = dose_character),
    size = 0.5,
    stroke = 0
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text =  element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_blank()
  ) +
  scale_color_manual(
    "Dose [nM]",
    values = c(
      "10" = "#1D1147FF",
      "100" = "#822681FF",
      "1000" = "#E65164FF",
      "10000" = "#FEC287FF"
    )
  ) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  ggsave(
    "supplemental_figure7_sciPlex_A549_24hr_UMAP_by_epoB_sub3.png",
    width = .75,
    height = .75,
    dpi = 900
  )


colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  filter(UMAP_1 > -6.5,
         UMAP_1 < -4.5,
         UMAP_2 < 1,
         UMAP_2 > -1) %>%
  group_by(product_name) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  head(n = 3)


colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  filter(UMAP_1 > 2,
         UMAP_1 < 4,
         UMAP_2 < -0.5,
         UMAP_2 > -2.5) %>%
  group_by(product_name) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  head(n = 3)

colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  filter(UMAP_1 > 3,
         UMAP_1 < 4,
         UMAP_2 < 2,
         UMAP_2 > 1) %>%
  group_by(product_name) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  head(n = 3)
##############


colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  group_by(Cluster) %>%
  add_tally() %>%
  ungroup() %>%
  dplyr::rename(cells_per_cluster = n) %>%
  group_by(catalog_number, Cluster) %>%
  mutate(percentage_cluster = n() / cells_per_cluster) %>%
  filter(catalog_number == "S1628") %>%
  dplyr::select(catalog_number, percentage_cluster, cells_per_cluster) %>%
  distinct()



colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  group_by(catalog_number) %>%
  add_tally() %>%
  ungroup() %>%
  dplyr::rename(cells_per_drug = n) %>%
  group_by(catalog_number, Cluster) %>%
  mutate(percentage_cluster = n() / cells_per_drug) %>%
  filter(catalog_number == "S1364") %>%
  dplyr::select(catalog_number, percentage_cluster, cells_per_drug) %>%
  distinct() %>%
  arrange(desc(percentage_cluster))

colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "grey80",
             size = 0.05,
             stroke = 0) +
  geom_point(
    data =
      colData(sciPlex_cds.list[["A549"]]) %>%
      as.data.frame() %>%
      filter(Cluster %in% c("9", "12", "16", "7")),
    aes(color = Cluster),
    size = 0.15,
    stroke = 0
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    text = element_text(size = 6),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.2, "line"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text =  element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_blank()
  ) +
  guides(guides(colour = guide_legend(override.aes = list(size = 1)))) +
  ggsave(
    "sciPlex_A549_clusters.png",
    width = .75,
    height = .75,
    dpi = 900
  )


colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  filter(Cluster %in% c("9", "12", "16", "7")) %>%
  group_by(pathway_level_1, Cluster) %>%
  add_tally() %>%
  top_n(3, wt = n) %>%
  dplyr::select(pathway_level_1, Cluster, n) %>%
  distinct() %>%
  arrange(desc(n))


colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  filter(Cluster %in% c("7")) %>%
  group_by(Cluster) %>%
  add_tally() %>%
  dplyr::rename(cells_in_cluster = n) %>%
  ungroup() %>%
  group_by(Cluster, product_name) %>%
  mutate(pathway_percent = n() / cells_in_cluster) %>%
  dplyr::select(cells_in_cluster, pathway_percent) %>%
  distinct() %>%
  arrange(desc(pathway_percent))

pathway_proportion_by_louvain_cluster.list <- list()

pathway_proportion_by_louvain_cluster.list[["A549"]] <-
  as.data.frame(
    pData(sciPlex_cds.list[["A549"]]) %>% as.data.frame() %>%
      filter(product_name != "Vehicle" , effective_dose == TRUE) %>%
      group_by(Cluster, pathway_level_1) %>%
      summarize(n = n()) %>%
      mutate(cluster_proportion = n /
               sum(n)) %>%
      dplyr::select(Cluster, pathway_level_1, cluster_proportion) %>%
      spread(key = pathway_level_1, value = cluster_proportion)
  )

pathway_proportion_by_louvain_cluster.list[["K562"]] <-
  as.data.frame(
    pData(sciPlex_cds.list[["K562"]]) %>% as.data.frame() %>%
      filter(product_name != "Vehicle" , effective_dose == TRUE) %>%
      group_by(Cluster, pathway_level_1) %>%
      summarize(n = n()) %>%
      mutate(cluster_proportion = n /
               sum(n)) %>%
      dplyr::select(Cluster, pathway_level_1, cluster_proportion) %>%
      spread(key = pathway_level_1, value = cluster_proportion)
  )

pathway_proportion_by_louvain_cluster.list[["MCF7"]] <-
  as.data.frame(
    pData(sciPlex_cds.list[["MCF7"]]) %>% as.data.frame() %>%
      filter(product_name != "Vehicle" , effective_dose == TRUE) %>%
      group_by(Cluster, pathway_level_1) %>%
      summarize(n = n()) %>%
      mutate(cluster_proportion = n /
               sum(n)) %>%
      dplyr::select(Cluster, pathway_level_1, cluster_proportion) %>%
      spread(key = pathway_level_1, value = cluster_proportion)
  )

cell_count_by_cluster.list <- list()

for (cell_line in names(pathway_proportion_by_louvain_cluster.list)) {
  row.names(pathway_proportion_by_louvain_cluster.list[[cell_line]]) <-
    pathway_proportion_by_louvain_cluster.list[[cell_line]]$Cluster
  
  pathway_proportion_by_louvain_cluster.list[[cell_line]]$Cluster <-
    NULL
  pathway_proportion_by_louvain_cluster.list[[cell_line]][is.na(pathway_proportion_by_louvain_cluster.list[[cell_line]])] <-
    0
  
  cell_count_by_cluster.list[[cell_line]] = sciPlex_cds.list[[cell_line]] %>%
    colData() %>%
    as.data.frame() %>%
    group_by(Cluster) %>%
    summarize(total_cells = log(n()))
  
}

pheatmap(
  pathway_proportion_by_louvain_cluster.list[["A549"]],
  file = "A549_drug_proportion_by_cluster_Residual_model_umi_replicate.png",
  width = 7.5,
  height = 10,
  annotation_row = data.frame(
    row.names = as.character(cell_count_by_cluster.list[["A549"]]$Cluster),
    total_cells = cell_count_by_cluster.list[["A549"]]$total_cells
  ),
  cluster_rows = F,
  color = plasma(35)
)

pheatmap(
  pathway_proportion_by_louvain_cluster.list[["K562"]],
  file = "K562_drug_proportion_by_cluster_Residual_model_umi_replicate.png",
  annotation_row = data.frame(
    row.names = as.character(cell_count_by_cluster.list[["K562"]]$Cluster),
    total_cells = cell_count_by_cluster.list[["K562"]]$total_cells
  ),
  width = 5.5,
  height = 10,
  cluster_rows = F,
  color = plasma(35)
)

pheatmap(
  pathway_proportion_by_louvain_cluster.list[["MCF7"]],
  file = "MCF7_drug_proportion_by_cluster_Residual_model_umi_replicate.png",
  annotation_row = data.frame(
    row.names = as.character(cell_count_by_cluster.list[["MCF7"]]$Cluster),
    total_cells = cell_count_by_cluster.list[["MCF7"]]$total_cells
  ),
  width = 5,
  height = 10,
  cluster_rows = F,
  color = plasma(35)
)




plot_percent_cells_positive(
  sciPlex_cds.list[["A549"]][fData(sciPlex_cds.list[["A549"]])$gene_short_name %in% c("GDF15", "ANGPTL4"),
                             pData(sciPlex_cds.list[["A549"]])$product_name %in% c("Vehicle",
                                                                                   "Triamcinolone Acetonide")],
  group_cells_by = "dose_character",
  plot_limits = NULL,
  ncol = 3
) +
  scale_fill_manual(
    "Log(Dose [µM])",
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
    axis.text.x = element_text(
      size = 4,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 4)
  ) +
  xlab("Log(Dose [µM])") +
  ylab("Cells (percent)") +
  ggsave(
    "supplemental_figure7_GR_marker_levels_in_Triamcinolone_exposed_A549_cells.png",
    width = 2,
    height = 0.8,
    unit = "in",
    dpi = 600
  )

plot_percent_cells_positive(
  sciPlex_cds.list[["K562"]][fData(sciPlex_cds.list[["K562"]])$gene_short_name %in% c("GDF15", "ANGPTL4"),
                             pData(sciPlex_cds.list[["K562"]])$product_name %in% c("Vehicle",
                                                                                   "Triamcinolone Acetonide")],
  group_cells_by = "dose_character",
  plot_limits = NULL,
  ncol = 3
) +
  scale_fill_manual(
    "Log(Dose [µM])",
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
    axis.text.x = element_text(
      size = 4,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 4)
  ) +
  xlab("Log(Dose [µM])") +
  ylab("Cells (percent)") +
  ggsave(
    "GR_marker_levels_in_Triamcinolone_exposed_K562cells.png",
    width = 2,
    height = 0.8,
    unit = "in",
    dpi = 600
  )

plot_percent_cells_positive(
  sciPlex_cds.list[["MCF7"]][fData(sciPlex_cds.list[["MCF7"]])$gene_short_name %in% c("GDF15", "ANGPTL4"),
                             pData(sciPlex_cds.list[["MCF7"]])$product_name %in% c("Vehicle",
                                                                                   "Triamcinolone Acetonide")],
  group_cells_by = "dose_character",
  plot_limits = NULL,
  ncol = 3
) +
  scale_fill_manual(
    "Log(Dose [µM])",
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
    axis.text.x = element_text(
      size = 4,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 4)
  ) +
  xlab("Log(Dose [µM])") +
  ylab("Cells (percent)") +
  ggsave(
    "GR_marker_levels_in_Triamcinolone_exposed_MCF7_cells.png",
    width = 2,
    height = 0.8,
    unit = "in",
    dpi = 600
  )
