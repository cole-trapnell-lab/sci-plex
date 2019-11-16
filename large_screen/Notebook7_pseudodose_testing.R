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
  library(tidyr)
  library(tidymodels)
  library(purrr)
  library(furrr)
  library(UpSetR)
  library(tictoc)
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
                "GSA_helper_functions.R",
                sep = ""))
  source(paste0(bin_directory,
                "loadGSCSafe.R",
                sep = ""))
})
# Set DelayedArray Parameters
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size = 1000e7)

# User defined functions

setwd(paste0(path_to_github,
             "large_screen",
             sep = ""))

cds = readRDS(file = "hdac_aligned_cds.RDS")

test_drug_dependence =
  function(cds_subset,
           cds,
           reference_cells = NULL,
           min_fraction_cells = 0.05,
           pseudodose = 0.01,
           residualModelFormulaTerms = NULL,
           cores = 1,
           return_model = F) {
    # For the subset of cells to get get their cell_ids
    cell_ids_to_test = as.character(cds_subset$cell)
    
    message(paste("Performing test on", length(cell_ids_to_test), "cells."))
    
    # Create a cds onject for the cell_ids that are going to be tested
    cds_subset = cds[, cell_ids_to_test]
    
    # Create a cds object for reference cells
    ref_cds = cds[, reference_cells]
    
    # Get the genes for the cells in the reference cds
    ref_cds = detect_genes(ref_cds)
    genes_in_reference_set = row.names(subset(
      fData(ref_cds),
      num_cells_expressed > ncol(ref_cds) * min_fraction_cells
    ))
    
    
    # Get the genes for the cells in the treatment subset cds
    cds_subset = detect_genes(cds_subset)
    genes_in_treated_set = row.names(subset(
      fData(cds_subset),
      num_cells_expressed > ncol(cds_subset) * min_fraction_cells
    ))
    
    message(paste(
      "Detected ",
      length(genes_in_treated_set),
      "genes for testing"
    ))
    
    cds_subset = cds[, base::union(cell_ids_to_test, row.names(pData(ref_cds)))]
    
    cds_subset = cds_subset[base::union(genes_in_reference_set, genes_in_treated_set), ]
    
    # At this point the cds_subset consists of the vehicle cells and the genes in each detected subset
    
    all_treatments =
      colData(cds_subset) %>%
      as.data.frame() %>%
      filter(!vehicle) %>%
      pull(treatment) %>%
      unique()
    
    for (tmnt in all_treatments) {
      cell_treated = (as.character(pData(cds_subset)$treatment) == as.character(tmnt))
      dose_of_drug = log10((pData(cds_subset)$dose * cell_treated)  + 1)
      drug_var_name = paste("dose_", tmnt, sep = "")
      pData(cds_subset)[, drug_var_name] = dose_of_drug
    }
    
    ct = colData(cds_subset)  %>%
      as.data.frame() %>%
      pull(cell_type) %>%
      unique()
    
    for (ct_type in ct) {
      print (ct_type)
      cell_treated = (as.character(pData(cds_subset)$cell_type) == as.character(ct_type))
      pData(cds_subset)[, ct_type] = as.numeric(cell_treated)
      message(paste(sum(cell_treated), ct_type, "in the dataset"))
    }
    
    message(paste("Fitting drug/dose models to ", nrow(cds_subset), "genes"))
    modelFormulaStr = "~0"
    
    
    for (ct_type in ct) {
      term = as.character(ct_type)
      #modelFormulaStr = paste(modelFormulaStr, "+" , term)
      modelFormulaStr = paste(modelFormulaStr,
                              " + " ,
                              term,
                              ":",
                              "splines::ns(Pseudotime, df=3)",
                              sep = "")
    }
    
    
    for (tmnt in all_treatments) {
      term = paste("dose_", tmnt, sep = "")
      modelFormulaStr = paste(modelFormulaStr, "+" , term)
    }
    
    print (modelFormulaStr)
    sci_plex_model_tbl = fit_models(
      cds_subset,
      model_formula_str = modelFormulaStr,
      cores = cores,
      reltol = 1e-5,
      verbose = TRUE
    )
    if (return_model == T) {
      return(sci_plex_model_tbl)
    } else{
      sci_plex_model_coefs = coefficient_table(sci_plex_model_tbl)
      return(sci_plex_model_coefs)
    }
    
  }

test_drug_dependence_model <-
  function(cds_subset,
           cds,
           reference_cells = NULL,
           min_fraction_cells = 0.05,
           pseudodose = 0.01,
           residualModelFormulaTerms = NULL,
           cores = 1) {
    # For the subset of cells to get get their cell_ids
    cell_ids_to_test = as.character(cds_subset$cell)
    
    message(paste("Performing test on", length(cell_ids_to_test), "cells."))
    
    # Create a cds onject for the cell_ids that are going to be tested
    cds_subset = cds[, cell_ids_to_test]
    
    # Create a cds object for reference cells
    ref_cds = cds[, reference_cells]
    
    # Get the genes for the cells in the reference cds
    ref_cds = detect_genes(ref_cds)
    genes_in_reference_set = row.names(subset(
      fData(ref_cds),
      num_cells_expressed > ncol(ref_cds) * min_fraction_cells
    ))
    
    
    # Get the genes for the cells in the treatment subset cds
    cds_subset = detect_genes(cds_subset)
    genes_in_treated_set = row.names(subset(
      fData(cds_subset),
      num_cells_expressed > ncol(cds_subset) * min_fraction_cells
    ))
    
    message(paste(
      "Detected ",
      length(genes_in_treated_set),
      "genes for testing"
    ))
    
    cds_subset = cds[, base::union(cell_ids_to_test, row.names(pData(ref_cds)))]
    #cds_subset = estimateDispersions(cds_subset)
    
    cds_subset = cds_subset[base::union(genes_in_reference_set, genes_in_treated_set), ]
    
    # At this point the cds_subset consists of the vehicle cells and the genes in each detected subset
    all_treatments = pData(cds_subset) %>%
      as.data.frame() %>%
      filter(!vehicle) %>%
      pull(treatment) %>%
      unique()
    
    for (tmnt in all_treatments) {
      #print (tmnt)
      cell_treated = (as.character(pData(cds_subset)$treatment) == as.character(tmnt))
      dose_of_drug = log10((pData(cds_subset)$dose * cell_treated)  + 1)
      drug_var_name = paste("dose_", tmnt, sep = "")
      pData(cds_subset)[, drug_var_name] = dose_of_drug
    }
    
    ct = pData(cds_subset)  %>%
      as.data.frame() %>%
      pull(cell_type) %>%
      unique()
    
    for (ct_type in ct) {
      print (ct_type)
      cell_treated = (as.character(pData(cds_subset)$cell_type) == as.character(ct_type))
      pData(cds_subset)[, ct_type] = as.numeric(cell_treated)
      message(paste(sum(cell_treated), ct_type, "in the dataset"))
    }
    
    message(paste("Fitting drug/dose models to ", nrow(cds_subset), "genes"))
    modelFormulaStr = "~0 +"
    
    
    for (ct_type in ct) {
      term = as.character(ct_type)
      modelFormulaStr = paste(modelFormulaStr, "+" , term)
      modelFormulaStr = paste(modelFormulaStr,
                              " + " ,
                              term,
                              ":",
                              "splines::ns(Pseudotime, df=3)",
                              sep = "")
    }
    
    
    for (tmnt in all_treatments) {
      term = paste("dose_", tmnt, sep = "")
      modelFormulaStr = paste(modelFormulaStr, "+" , term)
    }
    
    print (modelFormulaStr)
    
    sci_plex_model_tbl = fit_models(
      cds_subset,
      model_formula_str = modelFormulaStr,
      cores = cores,
      reltol = 1e-5,
      verbose = TRUE
    )
    return(sci_plex_model_tbl)
    
  }


cds = detect_genes(cds = cds, min_expr = 0.1)

# Method for picking tested_genes
genes_expressed_per_drug =
  colData(cds) %>%
  as.data.frame() %>%
  group_by(cell_type, treatment, time_point) %>%
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
    cds
  ))

genes_expressed_per_drug =
  genes_expressed_per_drug %>%
  unnest(fraction_genes)

expressed_genes =
  genes_expressed_per_drug %>%
  filter(fraction_cells_expressed > 0.01) %>%
  ungroup() %>%
  dplyr::select(id) %>%
  distinct() %>% pull(id)

cds_to_test = cds[expressed_genes, is.finite(colData(cds)$Pseudotime)]

vehicle_cells =
  colData(cds_to_test) %>%
  as.data.frame() %>%
  filter(vehicle == TRUE)

vehicle_cells = as.character(vehicle_cells$cell)


reference_cells =
  colData(cds_to_test) %>%
  as.data.frame() %>%
  sample_n(1000) %>%
  pull(cell)



pseudotime_model_res =
  colData(cds_to_test[expressed_genes, ]) %>%
  as.data.frame() %>%
  group_by() %>%
  tibble::rownames_to_column() %>%
  nest() %>%
  mutate(
    model_res = purrr::map(
      data,
      .f = test_drug_dependence_model,
      cds_to_test[expressed_genes, ],
      union(reference_cells, vehicle_cells),
      residualModelFormulaTerms = NULL,
      cores = 1
    )
  ) %>%
  dplyr::select(-data) %>%
  unnest()


pseudotime_coeff_res =
  colData(cds_to_test[expressed_genes, ]) %>%
  as.data.frame() %>%
  group_by() %>%
  tibble::rownames_to_column() %>%
  nest() %>%
  mutate(
    test_res_table = purrr::map(
      data,
      .f = test_drug_dependence,
      cds_to_test[expressed_genes, ],
      union(reference_cells, vehicle_cells),
      residualModelFormulaTerms = NULL,
      cores = 1
    )
  ) %>%
  dplyr::select(-data) %>%
  unnest()


model_terms =
  pseudotime_coeff_res %>%
  pull(term) %>%
  unique()

all_terms = c(
  "A549",
  "K562",
  "MCF7",
  "dose_S2693",
  "dose_S1045",
  "dose_S1096",
  "dose_S1053",
  "dose_S1030",
  "dose_S2170",
  "dose_S2759",
  "dose_S1085",
  "dose_S1095",
  "dose_S1194",
  "dose_S2818",
  "dose_S2244",
  "dose_S1090",
  "dose_S1515",
  "dose_S1122",
  "dose_S2779",
  "dose_S8567",
  "A549:splines::ns(Pseudotime, df = 3)1",
  "A549:splines::ns(Pseudotime, df = 3)2",
  "A549:splines::ns(Pseudotime, df = 3)3",
  "splines::ns(Pseudotime, df = 3)1:K562",
  "splines::ns(Pseudotime, df = 3)2:K562",
  "splines::ns(Pseudotime, df = 3)3:K562",
  "splines::ns(Pseudotime, df = 3)1:MCF7",
  "splines::ns(Pseudotime, df = 3)2:MCF7",
  "splines::ns(Pseudotime, df = 3)3:MCF7"
)

cell_line = c(
  "A549",
  "K562",
  "MCF7",
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  NA,
  "A549",
  "A549",
  "A549",
  "K562",
  "K562",
  "K562",
  "MCF7",
  "MCF7",
  "MCF7"
)


term_df = data.frame(term = all_terms,
                     cell_type = cell_line)

term_df$term_type = c(
  "Basal Expression",
  "Basal Expression",
  "Basal Expression",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Drug Dose",
  "Pseudodose A549",
  "Pseudodose A549",
  "Pseudodose A549",
  "Pseudodose K562",
  "Pseudodose K562",
  "Pseudodose K562",
  "Pseudodose MCF7",
  "Pseudodose MCF7",
  "Pseudodose MCF7"
)

term_df$drug = substr(term_df$term, start = 1, stop = 5) == "dose_"

term_df$Pseudotime = grepl("Pseudotime", term_df$term)

term_df$term_pretty <-
  plyr::revalue((term_df$term),
                c(
                  "A549" = "A549",
                  "K562" = "K562",
                  "MCF7" = "MCF7",
                  "dose_S2693" = "Resminostat",
                  "dose_S1045" = "Trichostatin",
                  "dose_S1096" = "Quisinostat",
                  "dose_S1053" = "Entinostat",
                  "dose_S1030" = "Panobinostat",
                  "dose_S2170" = "Givinostat",
                  "dose_S2759" = "CUDC-907",
                  "dose_S1085" = "Belinostat",
                  "dose_S1095" = "Dacinostat",
                  "dose_S1194" = "CUDC-101",
                  "dose_S2818" = "Tacedinaline",
                  "dose_S2244" = "AR-42",
                  "dose_S1090" = "Abexinostat",
                  "dose_S1515" = "Pracinostat",
                  "dose_S1122" = "Mocetinostat",
                  "dose_S2779" = "M344",
                  "dose_S8567" = "Tucidinostat",
                  "A549:splines::ns(Pseudotime, df = 3)1" = "Pseudodose1 * A549",
                  "A549:splines::ns(Pseudotime, df = 3)2" = "Pseudodose2 * A549",
                  "A549:splines::ns(Pseudotime, df = 3)3" = "Pseudodose3 * A549",
                  "splines::ns(Pseudotime, df = 3)1:K562" = "Pseudodose1 * K562",
                  "splines::ns(Pseudotime, df = 3)2:K562" = "Pseudodose2 * K562",
                  "splines::ns(Pseudotime, df = 3)3:K562" = "Pseudodose3 * K562",
                  "splines::ns(Pseudotime, df = 3)1:MCF7" = "Pseudodose1 * MCF7",
                  "splines::ns(Pseudotime, df = 3)2:MCF7" = "Pseudodose2 * MCF7",
                  "splines::ns(Pseudotime, df = 3)3:MCF7" = "Pseudodose3 * MCF7"
                )
  )

res_table =
  left_join(pseudotime_coeff_res, term_df, by = "term")



res_table %>%
  filter(q_value < .05) %>%
  filter(term_type != "Basal Expression") %>%
  dplyr::select(gene_short_name, term) %>%
  distinct() %>%
  group_by(term) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_bar(aes(x = reorder(term, n), y = n),
           color = "black",
           fill = "grey80",
           stat = "identity") +
  monocle:::monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.title.y = element_blank()
  ) +
  ylab("Number of DEGs") +
  coord_flip() +
  ggsave(
    "supplemental_figure28_differentially_expressed_gene_barplot_term.pdf",
    height = 1.25,
    width = 1.75
  )


pseudotime_collapsed =
  res_table %>%
  filter(Pseudotime,
         !is.na(cell_type),
         grepl("Pseudotime", term),
         q_value < .05) %>%
  ungroup()

pseudotime_collapsed =
  pseudotime_collapsed %>%
  mutate(Pseudotime_cell = paste(pseudotime_collapsed$cell_type, "Pseudodose", sep = " "))

pseudotime_terms =
  pseudotime_collapsed %>%
  pull(Pseudotime_cell) %>%
  unique()

listInput <- list()

for (t in pseudotime_terms) {
  print(t)
  listInput[[t]] =
    pseudotime_collapsed %>%
    filter(Pseudotime_cell == t) %>%
    pull(gene_short_name) %>%
    unique() %>%
    as.character()
}


pdf(file = "supplemental_figure28_hdac_trajectory_DEG_upsetr.pdf",
    onefile = FALSE,
    height = 2.5)
upset(fromList(listInput), order.by = "freq", nsets = 7)
dev.off()


# Plot the number of cells positive

# G2/M transition:
plot_percent_cells_positive(cds[fData(cds)$gene_short_name %in% c("CDKN1A", "AURKA"), ],
                            group_cells_by = "Pseudotime_bin",
                            plot_limits = NULL) +
  
  scale_fill_manual(values = c(
    "A549" = "#438FCD",
    "MCF7" = "#DB3C6A",
    "K562" = "#D6C939"
  )) +
  ggplot2::theme(legend.position = "none",
                 text = ggplot2::element_text(size = 6)) +
  scale_x_continuous(breaks = seq(2, 10, 2)) +
  xlab("Pseudodose Bin") +
  ggsave(
    "supplemental_figure29_hdac_pseudotime_g2m_markers.png",
    height = 1.25,
    width = 2,
    units = "in",
    dpi = 1200
  )


# Make Pseudodose Heatmaps
genes_to_test =
  pseudotime_coeff_res %>%
  filter(grepl("Pseudo", term), q_value < .05) %>%
  pull(id) %>%
  unique()


pseudotime_model_res =
  pseudotime_model_res %>%
  filter(id %in% genes_to_test)

min_pseudotime = min(colData(cds_to_test)$Pseudotime)
max_pseudotime = max(colData(cds_to_test)$Pseudotime)

A549_new_data = data.frame(
  Size_Factor = 1,
  cell_type = "A549",
  K562 = 0,
  A549 = 1,
  MCF7 = 0,
  dose_S2693 = 0,
  dose_S1045 = 0,
  dose_S1096 = 0,
  dose_S1053 = 0,
  dose_S1030 = 0,
  dose_S2170 = 0,
  dose_S2759 = 0,
  dose_S1085 = 0,
  dose_S1095 = 0,
  dose_S1194 = 0,
  dose_S2818 = 0,
  dose_S2244 = 0,
  dose_S1090 = 0,
  dose_S1515 = 0,
  dose_S1122 = 0,
  dose_S2779 = 0,
  dose_S8567 = 0,
  Pseudotime = seq(min_pseudotime, max_pseudotime, length.out = 100)
)

A549_preds = model_predictions(pseudotime_model_res, new_data = A549_new_data)
rownames(A549_preds) =  pseudotime_model_res$id

K562_new_data = data.frame(
  Size_Factor = 1 ,
  cell_type = "K562",
  K562 = 1,
  A549 = 0,
  MCF7 = 0,
  dose_S2693 = 0,
  dose_S1045 = 0,
  dose_S1096 = 0,
  dose_S1053 = 0,
  dose_S1030 = 0,
  dose_S2170 = 0,
  dose_S2759 = 0,
  dose_S1085 = 0,
  dose_S1095 = 0,
  dose_S1194 = 0,
  dose_S2818 = 0,
  dose_S2244 = 0,
  dose_S1090 = 0,
  dose_S1515 = 0,
  dose_S1122 = 0,
  dose_S2779 = 0,
  dose_S8567 = 0,
  Pseudotime = seq(min_pseudotime, max_pseudotime, length.out = 100)
)

K562_preds = model_predictions(pseudotime_model_res, new_data = K562_new_data)
rownames(K562_preds) =  pseudotime_model_res$id

MCF7_new_data = data.frame(
  Size_Factor = 1,
  cell_type = "MCF7",
  K562 = 0,
  A549 = 0,
  MCF7 = 1,
  dose_S2693 = 0,
  dose_S1045 = 0,
  dose_S1096 = 0,
  dose_S1053 = 0,
  dose_S1030 = 0,
  dose_S2170 = 0,
  dose_S2759 = 0,
  dose_S1085 = 0,
  dose_S1095 = 0,
  dose_S1194 = 0,
  dose_S2818 = 0,
  dose_S2244 = 0,
  dose_S1090 = 0,
  dose_S1515 = 0,
  dose_S1122 = 0,
  dose_S2779 = 0,
  dose_S8567 = 0,
  Pseudotime = seq(min_pseudotime, max_pseudotime, length.out = 100)
)
MCF7_preds = model_predictions(pseudotime_model_res, new_data = MCF7_new_data)
rownames(MCF7_preds) =  pseudotime_model_res$id


row_center = function(m) {
  # Row-center the data.
  m = m[!apply(m, 1, sd) == 0, ]
  m = Matrix::t(scale(Matrix::t(m), center = TRUE))
  m = m[is.na(row.names(m)) == FALSE, ]
  m[is.nan(m)] = 0
  return(m)
}

A549_preds = row_center(A549_preds)

K562_preds = row_center(K562_preds)

MCF7_preds = row_center(MCF7_preds)


predicted_heatmap = cbind(A549_preds, K562_preds, MCF7_preds)
rownames(predicted_heatmap) = pseudotime_model_res$id




heatmap_matrix <- predicted_heatmap

row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix))) / 2)
row_dist[is.na(row_dist)] <- 1


bks <- seq(-3.1, 3.1, by = 0.1)
hmcols <- viridis::viridis(length(bks) - 1)

heatmap_matrix[heatmap_matrix > 3] = 3

heatmap_matrix[heatmap_matrix < -3] = -3

ph <- pheatmap(
  heatmap_matrix,
  useRaster = T,
  cluster_cols = FALSE,
  cluster_rows = T,
  show_rownames = F,
  show_colnames = F,
  clustering_distance_rows = row_dist,
  clustering_method = "ward.D2",
  cutree_rows = 5,
  silent = TRUE,
  filename = NA,
  breaks = bks,
  color = hmcols
)

annotation_col = data.frame(CellType = factor(c(
  rep.int(x = "A549", times = 100),
  rep.int(x = "K562", times = 100),
  rep.int(x = "MCF7", times = 100)
)))



annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 5)))

col_gaps_ind <- c(1:(3 - 1)) * 100

ann_colors = list(
  Cluster = c(
    "1" = "#66c2a5",
    "2" = "#fc8d62",
    "3" = "#8da0cb",
    "4" = "#e78ac3",
    "5" = "#a6d854"
  ),
  CellType = c(
    "A549" = "#438FCD",
    "K562" = "#D6C939",
    "MCF7" = "#DB3C6A"
  )
)

ph_res <- pheatmap(
  heatmap_matrix,
  #ph$tree_row$order
  useRaster = T,
  cluster_cols = FALSE,
  cluster_rows = T,
  show_rownames = F,
  #scale="row",
  clustering_distance_rows = row_dist,
  clustering_method = "ward.D2",
  cutree_rows = 5,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  gaps_col = col_gaps_ind,
  treeheight_row = 0,
  fontsize = 12,
  color = hmcols,
  silent = F,
  border_color = NA,
  annotation_colors = ann_colors,
  annotation_legend = T,
  legend = F,
  annotation_names_row = F,
  annotation_names_col = F,
  filename = "supplemental_figure28_pseudotime_heatmap.png"
)

genes_to_plot =
  c(
    "SLC13A3",
    "SLC2A3",
    "GLS",
    "IDH1",
    "ACSS2",
    "ACSL3",
    "SIRT2",
    "SLC25A16",
    "SLC25A1",
    "ACLY",
    "HADHB",
    "HADHA",
    "ACO1",
    "CS"
  )

genes_to_plot_ids =
  rowData(cds) %>%
  as.data.frame() %>%
  filter(gene_short_name %in% genes_to_plot) %>%
  pull(id)

heatmap_matrix_subset = heatmap_matrix[genes_to_plot_ids, ]

rownames_new =
  rowData(cds) %>%
  as.data.frame() %>%
  filter(id %in% genes_to_plot_ids)


row.names(rownames_new) = rownames_new$id

rownames_new = rownames_new[rownames(heatmap_matrix_subset), ]
rownames(heatmap_matrix_subset) = rownames_new$gene_short_name

ph_res <- pheatmap(
  heatmap_matrix_subset,
  #ph$tree_row$order
  useRaster = T,
  cluster_cols = F,
  cluster_rows = T,
  show_rownames = T,
  show_colnames = F,
  #scale="row",
  clustering_method = "ward.D2",
  #ward.D2
  #cutree_cols = 3,
  gaps_col = col_gaps_ind,
  treeheight_row = 0,
  fontsize = 16,
  color = hmcols,
  silent = F,
  border_color = NA,
  legend = F,
  annotation_names_row = F,
  annotation_names_col = F,
  filename = "figure5_pseudotime_heatmap_subset.png"
)

write.table(
  x = pseudotime_coeff_res,
  file = "Supplementary_Table_8_HDACi_DEGs.txt",
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)

######################




plot_percent_aligned_cells_positive <- function(cds_subset,
                                                grouping = "Pseudotime_bin",
                                                min_expr = 0.1,
                                                nrow = NULL,
                                                ncol = 1,
                                                panel_order = NULL,
                                                plot_as_fraction = TRUE,
                                                label_by_short_name = TRUE,
                                                relative_expr = TRUE,
                                                plot_limits = NULL) {
  percent <- NULL
  
  integer_expression <- TRUE
  
  if (integer_expression)
  {
    marker_exprs <- exprs(cds_subset)
    if (relative_expr) {
      if (is.null(size_factors(cds_subset)))
      {
        stop(
          "Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first"
        )
      }
      marker_exprs <-
        Matrix::t(Matrix::t(marker_exprs) / size_factors(cds_subset))
    }
    marker_exprs_melted <-
      reshape2::melt(round(as.matrix(marker_exprs)))
  } else{
    marker_exprs_melted <- reshape2::melt(exprs(marker_exprs))
  }
  
  colnames(marker_exprs_melted) <- c("f_id", "Cell", "expression")
  
  marker_exprs_melted <-
    merge(marker_exprs_melted,
          as.data.frame(pData(cds_subset)),
          by.x = "Cell",
          by.y = "row.names")
  marker_exprs_melted <-
    merge(marker_exprs_melted,
          as.data.frame(fData(cds_subset)),
          by.x = "f_id",
          by.y = "row.names")
  
  if (label_by_short_name == TRUE) {
    if (is.null(marker_exprs_melted$gene_short_name) == FALSE) {
      marker_exprs_melted$feature_label <-
        marker_exprs_melted$gene_short_name
      marker_exprs_melted$feature_label[is.na(marker_exprs_melted$feature_label)]  <-
        marker_exprs_melted$f_id
    } else{
      marker_exprs_melted$feature_label <- marker_exprs_melted$f_id
    }
  } else{
    marker_exprs_melted$feature_label <- marker_exprs_melted$f_id
  }
  
  if (is.null(panel_order) == FALSE)
  {
    marker_exprs_melted$feature_label <-
      factor(marker_exprs_melted$feature_label, levels = panel_order)
  }
  
  marker_counts <-
    plyr::ddply(marker_exprs_melted, c("feature_label", "cell_type", grouping), function(x) {
      data.frame(
        target = sum(x$expression > min_expr),
        target_fraction = sum(x$expression > min_expr) / nrow(x)
      )
    })
  
  #print (head(marker_counts))
  if (plot_as_fraction) {
    marker_counts$target_fraction <- marker_counts$target_fraction * 100
    qp <-
      ggplot2::ggplot(ggplot2::aes_string(x = grouping, y = "target_fraction", fill =
                                            "cell_type"),
                      data = marker_counts) +
      ggplot2::ylab("Cells (percent)")
    if (is.null(plot_limits) == FALSE)
      qp <- qp + ggplot2::scale_y_continuous(limits = plot_limits)
  } else{
    qp <-
      ggplot2::ggplot(ggplot2::aes_string(x = grouping, y = "target", fill =
                                            "cell_type"),
                      data = marker_counts) +
      ggplot2::ylab("Cells")
  }
  
  qp <-
    qp + ggplot2::facet_grid(feature_label ~ cell_type, scales = "free_y")
  qp <-
    qp + ggplot2::geom_bar(stat = "identity",
                           color = "black",
                           size = .2) + monocle3:::monocle_theme_opts()
  
  return(qp)
}



# G2/M transition:
plot_percent_cells_positive(cds[fData(cds)$gene_short_name %in% c("CDKN1A", "AURKA"), pData(cds)$cell_type ==
                                  "A549"],
                            group_cells_by = "Pseudotime_bin") +
  scale_fill_viridis() +
  ggplot2::theme(
    legend.position = "none",
    text = ggplot2::element_text(size = 6),
    axis.title = element_blank(),
    strip.text = element_blank(),
    axis.ticks = element_line(size = 0.15)
  ) +
  scale_x_continuous(breaks = seq(2, 10, 2)) +
  ggsave(
    "supplemental_figure29_a549_hdac_pseudotime_g2m_markers.pdf",
    height = 1,
    width = .75,
    units = "in"
  )


# G2/M transition:
plot_percent_cells_positive(cds[fData(cds)$gene_short_name %in% c("CDKN1A", "AURKA"), pData(cds)$cell_type ==
                                  "K562"],
                            group_cells_by = "Pseudotime_bin") +
  scale_fill_viridis() +
  ggplot2::theme(
    legend.position = "none",
    text = ggplot2::element_text(size = 6),
    axis.title = element_blank(),
    strip.text = element_blank(),
    axis.ticks = element_line(size = 0.15)
  ) +
  scale_x_continuous(breaks = seq(2, 10, 2)) +
  ggsave(
    "supplemental_figure29_k562_hdac_pseudotime_g2m_markers.pdf",
    height = 1,
    width = .75,
    units = "in"
  )

plot_percent_cells_positive(cds[fData(cds)$gene_short_name %in% c("CDKN1A", "AURKA"), pData(cds)$cell_type ==
                                  "MCF7"],
                            group_cells_by = "Pseudotime_bin") +
  scale_fill_viridis() +
  ggplot2::theme(
    legend.position = "none",
    text = ggplot2::element_text(size = 6),
    axis.title = element_blank(),
    strip.text = element_blank(),
    axis.ticks = element_line(size = 0.15)
  ) +
  scale_x_continuous(breaks = seq(2, 10, 2)) +
  ggsave(
    "supplemental_figure29_mcf7_hdac_pseudotime_g2m_markers.pdf",
    height = 1,
    width = .75,
    units = "in"
  )



#####################

## Load Gene Set Collection
hallmarks <- loadGSCSafe(file = paste0(path_to_github,
                                       "GeneSets/h.all.v6.0.symbols.gmt",
                                       sep = ""))

genes_for_GSAhyper = annotation_row
genes_for_GSAhyper$id  = rownames(genes_for_GSAhyper)

genes_for_GSAhyper = left_join(genes_for_GSAhyper, rowData(cds) %>% as.data.frame(), by = "id")


clust1 =
  genes_for_GSAhyper %>%
  filter(Cluster == "1") %>%
  pull(gene_short_name) %>%
  unique()

clust1 = unique(clust1)
gsaRes1 <-
  runGSAhyper(clust1, gsc = hallmarks, universe = gene_universe)


clust2 =
  genes_for_GSAhyper %>%
  filter(Cluster == "2") %>%
  pull(gene_short_name) %>%
  unique()


clust2 = unique(clust2)
gsaRes2 <-
  runGSAhyper(clust2, gsc = hallmarks, universe = gene_universe)


clust3 =
  genes_for_GSAhyper %>%
  filter(Cluster == "3") %>%
  pull(gene_short_name) %>%
  unique()


clust3 = unique(clust3)
gsaRes3 <-
  runGSAhyper(clust3, gsc = hallmarks, universe = gene_universe)

clust4 =
  genes_for_GSAhyper %>%
  filter(Cluster == "4") %>%
  pull(gene_short_name) %>%
  unique()


clust4 = unique(clust4)
gsaRes4 <-
  runGSAhyper(clust4, gsc = hallmarks, universe = gene_universe)


clust5 =
  genes_for_GSAhyper %>%
  filter(Cluster == "5") %>%
  pull(gene_short_name) %>%
  unique()


clust5 = unique(clust5)
gsaRes5 <-
  runGSAhyper(clust5, gsc = hallmarks, universe = gene_universe)


collect_gsa_hyper_results <- function(cds, gsc, clusters)
{
  gene_universe <- unique(as.character(fData(cds)$gene_short_name))
  gsa_results <- list()
  cluster_ids <- unique(clusters)
  for (i in (1:length(cluster_ids))) {
    cluster_genes <-
      unique(fData(cds[names(clusters[clusters == i]), ])$gene_short_name)
    gsaRes <-
      runGSAhyper(cluster_genes, gsc = gsc, universe = gene_universe)
    gsa_results[[length(gsa_results) + 1]] <- gsaRes
  }
  names(gsa_results) <- cluster_ids
  gsa_results
}