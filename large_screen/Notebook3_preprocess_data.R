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
  library(reshape2)
  library(Matrix)
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
})

# Set DelayedArray Parameters
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size = 1000e7)

# User defined functions

setwd(paste0(path_to_github,
             "large_screen",
             sep = ""))


# Read in CDS and extract cells passing filters and human transcripts
cds = readRDS(cds.RDS)

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


# Isolate cells at 24 hours for initial analysis
cells_24hrs =
  colData(cds) %>%
  as.data.frame() %>%
  filter(time_point == 24) %>%
  pull(cell)

# Subset 24 hour cells
cds_24hrs = cds[, cells_24hrs]

# Pull out the cells corresponding to each cell type
sciPlex_cds.list <- list()

for (this.cell_type in unique(colData(cds_24hrs)$cell_type)) {
  these.cells =
    colData(cds_24hrs) %>%
    as.data.frame() %>%
    filter(cell_type == this.cell_type) %>%
    pull(cell)
  
  this.cds = cds_24hrs[, these.cells]
  
  sciPlex_cds.list[[this.cell_type]] = this.cds
  
  print(paste0(
    "Recovered ",
    this.cell_type,
    " ",
    length(these.cells),
    " cells"
  ))
}

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
  
  # Fit dose response curves based on cell counts
  sciPlex_models.list[[cell_line]] =
    fit_dose_response_models(sciPlex_counts.list[[cell_line]])
  
  # Recover model statistics for the dose response curves
  sciPlex_model_stats.list[[cell_line]] =
    extract_dose_model_stats(sciPlex_counts.list[[cell_line]], sciPlex_models.list[[cell_line]])
  
  colData(sciPlex_cds.list[[cell_line]])$product_dose =
    paste0(colData(sciPlex_cds.list[[cell_line]])$product_name,
           "_",
           colData(sciPlex_cds.list[[cell_line]])$dose)
}


# Get genes that are expressed in greater than 5% of cells in a treatment or timepoint
# Subset of genes to be used for PCA and results that use PCA results
expressed_genes.list <- list()

for (cell_line in names(sciPlex_cds.list)) {
  genes_expressed_per_drug = colData(sciPlex_cds.list[[cell_line]]) %>%
    as.data.frame() %>%
    group_by(treatment, time_point) %>%
    nest() %>%
    mutate(fraction_genes = purrr::map(
      data,
      .f = function(pdata_subset, cds) {
        cds_subset = cds[, as.character(pdata_subset$cell)]
        cds_subset = detect_genes(cds_subset)
        tibble(
          id = rowData(cds_subset)$id,
          gene_short_name = rowData(cds_subset)$gene_short_name,
          num_cells_expressed = rowData(cds_subset)$num_cells_expressed,
          fraction_cells_expressed = rowData(cds_subset)$num_cells_expressed / ncol(cds_subset)
        )
      },
      sciPlex_cds.list[[cell_line]]
    ))
  
  genes_expressed_per_drug = genes_expressed_per_drug %>% unnest(fraction_genes)
  expressed_genes = genes_expressed_per_drug %>% filter(fraction_cells_expressed > 0.05) %>% ungroup() %>% dplyr::select(id) %>% distinct()
  expressed_genes.list[[cell_line]] = as.character(expressed_genes$id)
  
  print(length(expressed_genes.list[[cell_line]]))
}


### For each cell line, perform the following:
# 1) Viability estimation
# 2) Detect genes
# 3) Estimate size factors
# 4) Get top 25 PCs
# 5) Perform UMAP dimensionality reduction
# 6) Cluster Cells in PCA space amd in UMAP space
resolution_louvain = list()
resolution_louvain[["A549"]] = 0.0005
resolution_louvain[["MCF7"]] = 0.00025
resolution_louvain[["K562"]] = 0.00025
# 7) Learn a graph in partitioned UMAP space


for (cell_line in names(sciPlex_cds.list)) {
  # Perform viability estimation on each cell line with the fit sciPlex models
  sciPlex_cds.list[[cell_line]] =
    estimate_viability(sciPlex_cds.list[[cell_line]], sciPlex_models.list[[cell_line]])
  
  # Set all viability that are estimated as NA to 1
  colData(sciPlex_cds.list[[cell_line]])$viability[is.na(pData(sciPlex_cds.list[[cell_line]])$viability)] = 1
  
  # Detect the genes that are expressed in the dataset
  sciPlex_cds.list[[cell_line]] <-
    detect_genes(sciPlex_cds.list[[cell_line]], min_expr = 0.1)
  
  # Estimate Size Factors for every cell
  sciPlex_cds.list[[cell_line]] <-
    estimate_size_factors(sciPlex_cds.list[[cell_line]])
  
  # Estimate Cell Cycle based on marker gene expression for every cell
  sciPlex_cds.list[[cell_line]] <-
    estimate_cell_cycle(
      sciPlex_cds.list[[cell_line]],
      g1s_markers = cc.genes$s.genes,
      g2m_markers = cc.genes$g2m.genes
    )
  
  
  # After regresing out the effects of number of UMIs and replicate, perform PCA and return 25 PCs
  sciPlex_cds.list[[cell_line]] <-
    preprocess_cds(
      sciPlex_cds.list[[cell_line]],
      method = 'PCA',
      num_dim = 25,
      norm_method = 'log',
      use_genes = expressed_genes.list[[cell_line]],
      residual_model_formula_str = "~log(n.umi) + replicate",
      verbose = T,
      cores = 20
    )
  
  fname = paste(cell_line, "_24hrs_scree.pdf", sep = "")
  monocle3::plot_pc_variance_explained(cds = sciPlex_cds.list[[cell_line]]) +
    ggsave(fname)
  
  # Perform UMAP dimensionality reduction
  sciPlex_cds.list[[cell_line]] <-
    reduce_dimension(
      sciPlex_cds.list[[cell_line]],
      max_components = 2,
      reduction_method = 'UMAP',
      umap.metric = 'cosine',
      umap.n_neighbors = 20,
      umap.min_dist = 0.1,
      verbose = TRUE,
      cores = 20
    )
  
  # Append UMAP coordinates to the CDS colData
  colData(sciPlex_cds.list[[cell_line]])$UMAP_1 <-
    sciPlex_cds.list[[cell_line]]@reducedDims[["UMAP"]][, 1]
  colData(sciPlex_cds.list[[cell_line]])$UMAP_2 <-
    sciPlex_cds.list[[cell_line]]@reducedDims[["UMAP"]][, 2]
  
  
  # Cluster based on the PCA dimensionality reduced data
  sciPlex_cds.list[[cell_line]] <-
    cluster_cells(sciPlex_cds.list[[cell_line]],
                  reduction_method = "PCA",
                  resolution = resolution_louvain[[cell_line]])
  
  # Append cluster information to colData
  colData(sciPlex_cds.list[[cell_line]])$Cluster = clusters(sciPlex_cds.list[[cell_line]], reduction_method =
                                                              "PCA")
  
  sciPlex_cds.list[[cell_line]] <-
    cluster_cells(sciPlex_cds.list[[cell_line]],
                  reduction_method = "UMAP")
  
  colData(sciPlex_cds.list[[cell_line]])$louvain_component = sciPlex_cds.list[[cell_line]]@clusters[["UMAP"]]$partitions
  
  sciPlex_cds.list[[cell_line]] = learn_graph(sciPlex_cds.list[[cell_line]])
  
  file_name = paste(cell_line, "_24hrs", ".RDS", sep = "")
  saveRDS(object = sciPlex_cds.list[[cell_line]], file =  file_name)
  
}
