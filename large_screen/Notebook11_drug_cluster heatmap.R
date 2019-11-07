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
  library(tidymodels)
  library(devtools)
  load_all(path_to_monocle3)
})



# Set DelayedArray Parameters
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size = 1000e7)

# User defined functions

setwd(paste0(path_to_github,
             "large_screen",
             sep = ""))

product_cluster_mat.list <- list()
sciPlex_cds.list = list(
  "A549" = readRDS("A549_24hrs.RDS"),
  "K562" = readRDS("K562_24hrs.RDS"),
  "MCF7" = readRDS("MCF7_24hrs.RDS")
)

###########

for (cell_line in names(sciPlex_cds.list)) {
  product_cluster_mat.list[[cell_line]] = reshape2::acast(
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
}


cluster.enrichment.df = list()

for (cell_line in names(sciPlex_cds.list)) {
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
}


for (cell_line in names(sciPlex_cds.list)) {
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
#########

min_enrichment_over_vehicle = -2
max_enrichment_over_vehicle = 2

drug_heatmap_annotations = bind_rows(
  colData(sciPlex_cds.list[["A549"]]) %>% as.data.frame,
  colData(sciPlex_cds.list[["K562"]]) %>% as.data.frame,
  colData(sciPlex_cds.list[["MCF7"]]) %>% as.data.frame
) %>%
  select(pathway_level_2, pathway_level_1, product_name) %>% distinct() %>% arrange(pathway_level_1, pathway_level_2)
drug_labels = drug_heatmap_annotations %>% pull(product_name)
row.names(drug_heatmap_annotations) = drug_heatmap_annotations$product_name
drug_heatmap_annotations$product_name = NULL

cluster_labels = c()
for (cell_line in names(sciPlex_cds.list)) {
  cluster_labels = c(cluster_labels, paste(cell_line, unique(colData(
    sciPlex_cds.list[[cell_line]]
  )$Cluster)))
}

cluster_enrich_over_vehicle = matrix(NA,
                                     nrow = length(drug_labels),
                                     ncol = length(cluster_labels))
colnames(cluster_enrich_over_vehicle) = cluster_labels
row.names(cluster_enrich_over_vehicle) = drug_labels

cluster_heatmap_annotations = data.frame()

for (cell_line in names(sciPlex_cds.list)) {
  print(cell_line)
  cds = sciPlex_cds.list[[cell_line]]
  
  vehicle_cell_ids = colData(cds) %>% as.data.frame %>% rownames_to_column %>% filter(product_name == "Vehicle") %>% pull(cell)
  vehicle_cell_ids = as.character(vehicle_cell_ids)
  vehicle_cds = cds[, vehicle_cell_ids]
  
  vehicle_freqs = as.numeric(prop.table(table(colData(cds[, vehicle_cell_ids])$Cluster)))
  drug_freqs = as.matrix(prop.table(table(
    colData(cds)$product_name, colData(cds)$Cluster
  ), margin = 1))
  
  cell_line_cluster_enrich_over_vehicle = t(t(drug_freqs) / vehicle_freqs)
  cell_line_cluster_enrich_over_vehicle = log2(cell_line_cluster_enrich_over_vehicle)
  cell_line_cluster_enrich_over_vehicle[cell_line_cluster_enrich_over_vehicle < min_enrichment_over_vehicle] = min_enrichment_over_vehicle
  cell_line_cluster_enrich_over_vehicle[cell_line_cluster_enrich_over_vehicle > max_enrichment_over_vehicle] = max_enrichment_over_vehicle
  cell_line_cluster_enrich_over_vehicle[is.na(cell_line_cluster_enrich_over_vehicle)] = 0
  
  cluster_drug_combos = cluster.enrichment.df[[cell_line]] %>% group_by(product_name, Cluster) %>% filter(q.value < 0.05) %>% top_n(1,-q.value) %>% select(product_name, Cluster) %>% distinct()
  cluster_drug_combos = cluster_drug_combos %>%
    mutate(yesno = 1) %>%
    distinct %>%
    spread(Cluster, yesno, fill = NA)
  cluster_drug_combos = as.data.frame(cluster_drug_combos)
  row.names(cluster_drug_combos) = cluster_drug_combos$product_name
  cluster_drug_combos$product_name = NULL
  cluster_drug_combos = as.matrix(cluster_drug_combos)
  cell_line_cluster_enrich_over_vehicle = cell_line_cluster_enrich_over_vehicle[row.names(cluster_drug_combos), colnames(cluster_drug_combos)]
  cell_line_cluster_enrich_over_vehicle = cell_line_cluster_enrich_over_vehicle * cluster_drug_combos
  
  colnames(cell_line_cluster_enrich_over_vehicle) = paste(cell_line,
                                                          colnames(cell_line_cluster_enrich_over_vehicle))
  cluster_enrich_over_vehicle[row.names(cell_line_cluster_enrich_over_vehicle), colnames(cell_line_cluster_enrich_over_vehicle)] = cell_line_cluster_enrich_over_vehicle
  
  cell_line_cluster_heatmap_annotations = colData(cds) %>% as.data.frame %>% group_by(Cluster) %>%
    summarize(
      proliferation = median(proliferation_index),
      viability = mean(viability)
    ) %>% as.data.frame
  cell_line_cluster_heatmap_annotations$Cluster = as.character(cell_line_cluster_heatmap_annotations$Cluster)
  cell_line_cluster_heatmap_annotations$Cluster = paste(cell_line, cell_line_cluster_heatmap_annotations$Cluster)
  cluster_heatmap_annotations = bind_rows(cluster_heatmap_annotations,
                                          cell_line_cluster_heatmap_annotations)
}

cluster_heatmap_annotations$cell_line = stringr::str_split_fixed(cluster_heatmap_annotations$Cluster, " ", 2)[, 1]
row.names(cluster_heatmap_annotations) = cluster_heatmap_annotations$Cluster
cluster_heatmap_annotations$Cluster = NULL

breaksList = seq(min_enrichment_over_vehicle,
                 max_enrichment_over_vehicle,
                 by = 0.1)

drug_order = colData(cds) %>% as.data.frame %>% select(pathway_level_1, pathway_level_2, target, product_name) %>% distinct() %>% arrange(pathway_level_1, pathway_level_2, target) %>% pull(product_name)

cluster_drugs = function(product_names, enrichment_mat) {
  na_mask = is.na(enrichment_mat)
  enrichment_mat[na_mask] = 0
  print (product_names)
  if (length(product_names) > 2) {
    cluster_res = hclust(dist(enrichment_mat[product_names, ]), method = "ward.D2")
    data.frame(drug_order = as.character(cluster_res$labels[cluster_res$order]))
  } else{
    data.frame(drug_order = product_names)
  }
}
drug_order = colData(cds) %>% as.data.frame %>% select(pathway_level_1, pathway_level_2, product_name) %>%
  distinct() %>% group_by(pathway_level_1, pathway_level_2) %>% do(cluster_drugs(.$product_name, cluster_enrich_over_vehicle))
drug_order = drug_order %>% pull(drug_order)


cluster_responses = function(response_names, enrichment_mat) {
  na_mask = is.na(enrichment_mat)
  enrichment_mat[na_mask] = 0
  print (response_names)
  if (length(response_names) > 2) {
    cluster_res = hclust(dist(t(enrichment_mat[, response_names])), method =
                           "ward.D2")
    data.frame(response_order = as.character(cluster_res$labels[cluster_res$order]))
  } else{
    data.frame(response_order = product_names)
  }
}
response_order = cluster_heatmap_annotations %>% rownames_to_column() %>% group_by(cell_line) %>% do(cluster_responses(.$rowname, cluster_enrich_over_vehicle))
response_order = response_order %>% pull(response_order)


na_mask = is.na(cluster_enrich_over_vehicle)
cluster_enrich_over_vehicle[na_mask] = 0
col_clust = hclust(dist(t(cluster_enrich_over_vehicle)))
cluster_enrich_over_vehicle[na_mask] = NA

cluster_enrich_over_vehicle = cluster_enrich_over_vehicle[drug_order, response_order]
drug_heatmap_annotations = drug_heatmap_annotations[drug_order, ]
cluster_heatmap_annotations = cluster_heatmap_annotations[response_order, ]

drug_heatmap_annotations$pathway_level_2 = stringr::str_c(
  drug_heatmap_annotations$pathway_level_1,
  " - ",
  drug_heatmap_annotations$pathway_level_2
)


column_gaps = cluster_heatmap_annotations %>% mutate(row_idx = row_number()) %>% group_by(cell_line) %>% top_n(1) %>% pull(row_idx)
row_gaps = drug_heatmap_annotations %>% mutate(row_idx = row_number()) %>% group_by(pathway_level_1) %>% top_n(1) %>% pull(row_idx)


pheatmap::pheatmap(
  cluster_enrich_over_vehicle,
  fontsize = 4,
  cluster_rows = FALSE,
  #cluster_cols = col_clust,
  cluster_cols = FALSE,
  gaps_col = column_gaps,
  gaps_row = row_gaps,
  annotation_row = drug_heatmap_annotations,
  annotation_col = cluster_heatmap_annotations,
  filename = "supplemental_figure6_cluster_enrich_vs_vehicle.pdf",
  color = colorRampPalette(rev(
    RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")
  ))(length(breaksList)),
  # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
  breaks = breaksList,
  height = 10,
  width = 6
)
