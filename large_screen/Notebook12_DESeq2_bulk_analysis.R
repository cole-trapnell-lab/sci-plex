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
  library(stringr)
  library(tidyr)
  library(purrr)
  library(ggridges)
  library(tidymodels)
  library(snowfall)
  library(DESeq2)
  library(ggpubr)
  library(devtools)
  library(piano)
  load_all(path_to_monocle3)
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

coldata_cds =
  colData(cds) %>%
  as.data.frame() %>%
  mutate(
    rt_index = stringr::str_split_fixed(Cell, "RT_", 2)[, 2] %>% as.character(),
    PCR_well = stringr::str_sub(string = Cell, 1, 1),
    cleanup = ifelse(PCR_well %in% c("A", "B", "C", "D"), "SPRI", "Exo")
  )


count_mat_merged =
  coldata_cds %>%
  group_by(rt_index) %>%
  nest() %>%
  mutate(counts = purrr::map(
    data,
    .f = function(data_subset, cds) {
      cds_subset = cds[, data_subset$Cell]
      totals = Matrix::rowSums(counts(cds_subset))
      tibble(gene_id = rownames(exprs(cds_subset)),
             counts = totals)
    },
    cds
  )) %>%
  dplyr::select(rt_index, counts) %>%
  unnest() %>%
  mutate(rt_index = paste("RT_", rt_index, sep = "")) %>%
  spread(key = rt_index, value = counts)

count_mat_merged =
  count_mat_merged %>%
  tibble::column_to_rownames(var = "gene_id") %>%
  as.matrix()

coldata = data.frame(row.names = colnames(count_mat_merged),
                     cell = colnames(count_mat_merged))

coldata_rts = read.table("bin/drug_rt_indices.tsv", header = F)

colnames(coldata_rts) = c("rt_index", "catalog_number", "dose", "cell_type")

coldata_rts$rt_index =
  paste("RT_", coldata_rts$rt_index, sep = "") %>%
  as.character()

coldata$cell = as.character(coldata$cell)

coldata = left_join(coldata, coldata_rts, by = c("cell" = "rt_index"))
rownames(coldata) = coldata$cell


coldata$catalog_number = factor(coldata$catalog_number)
coldata$catalog_number = relevel(coldata$catalog_number, ref  = "S0000")


rowdata = data.frame(
  row.names = rownames(count_mat_merged),
  id = rownames(count_mat_merged) %>%
    as.character()
)


############ All gene TPM correlations  #############

coldata$drug_dose =
  paste(
    coldata$cell_type %>% as.character(),
    coldata$catalog_number %>% as.character(),
    coldata$dose %>% as.character(),
    "bulk",
    sep = "_"
  ) %>%
  as.character()

coldata$library = "bulk"

# Path to cds containing all the data
path_to_cds = "cds.RDS"

# Read in CDS and extract cells passing filters and human transcripts
cds_sciPlex = readRDS(path_to_cds)

# Genes in the human transcriptome
human.genes = grepl("^ENSG", rowData(cds_sciPlex)$id)


# Cells that passed hash filters
validly_labeled_cells =
  colData(cds_sciPlex) %>%
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

cds_sciPlex = cds_sciPlex[human.genes, validly_labeled_cells]

# Isolate cells at treated with subset of drugs for 24 hours for analysis
cells_to_pull =
  colData(cds_sciPlex) %>%
  as.data.frame() %>%
  filter(
    time_point == 24,
    cell_type %in% c("K562", "A549"),
    catalog_number %in% coldata$catalog_number
  ) %>%
  pull(cell)

cds_sciPlex = cds_sciPlex[, cells_to_pull]

# Merge_cds
coldata_sci =
  colData(cds_sciPlex) %>%
  as.data.frame() %>%
  dplyr::select("cell", "catalog_number", "dose", "cell_type") %>%
  mutate(drug_dose = paste(cell_type,
                           catalog_number,
                           dose %>% as.character(),
                           "sci",
                           sep = "_")) %>%
  mutate(library = "sci")


coldata_bulk = coldata

row_order = rownames(cds_sciPlex)

count_mat_merged = count_mat_merged[row_order, ]

joint_matrix = cbind(as.matrix(counts(cds_sciPlex)), count_mat_merged)
joint_coldata = rbind(coldata_sci, coldata_bulk)
rownames(joint_coldata) = joint_coldata$cell

joint_rowdata = data.frame(row.names = rownames(joint_matrix),
                           id = rownames(joint_matrix))

joint_cds = new_cell_data_set(joint_matrix,
                              joint_coldata,
                              joint_rowdata)

joint_cds = estimate_size_factors(joint_cds)

joint_matrix_aggregated =
  colData(joint_cds) %>%
  as.data.frame() %>%
  group_by(drug_dose) %>%
  nest() %>%
  mutate(aggregated_mat = purrr::map(data, function(coldata_subset, cds) {
    cell_subset = coldata_subset$cell
    count_matrix_subset = cds[, cell_subset]
    norm_exprs_mat = Matrix::t(Matrix::t(counts(count_matrix_subset)) /
                                 size_factors(count_matrix_subset))
    totals = Matrix::rowSums(counts(count_matrix_subset))
    tibble(gene_id = rownames(count_matrix_subset),
           counts = totals)
  }, joint_cds)) %>%
  dplyr::select(-data) %>%
  unnest() %>%
  ungroup()

joint_matrix_aggregated =
  joint_matrix_aggregated %>%
  tidyr::spread(key = drug_dose, value = counts,-2) %>%
  tibble::column_to_rownames(var = "gene_id") %>%
  as.matrix()

joint_matrix_aggregated = Matrix::t(Matrix::t(joint_matrix_aggregated) /
                                      Matrix::colSums(joint_matrix_aggregated)) * 1000000 + 1


joint_matrix_aggregated =
  joint_matrix_aggregated %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "id") %>%
  gather(key = drug_dose, value = counts,-id)


joint_matrix_aggregated =
  joint_matrix_aggregated %>%
  separate(
    col = drug_dose,
    into = c("cell_type", "catalog_number", "dose", "library"),
    sep = "_"
  ) %>%
  mutate(unique_id = paste(cell_type, catalog_number, dose, sep = "_")) %>%
  dplyr::select(unique_id, id, library, counts) %>%
  spread(key = library, value = counts,-1:2)

joint_matrix_aggregated =
  joint_matrix_aggregated %>%
  separate(unique_id,
           into = c("cell_type", "drug", "dose"),
           sep = "_")

rowdata_genes =
  read.table(
    file = "bin/gencode.v27.transcripts.bed",
    header = F,
    sep = "\t",
    col.names = c(
      "chromosome",
      "start",
      "end",
      "tag",
      "dummy",
      "strand",
      "id",
      "gene_short_name",
      "classification"
    )
  ) %>%
  dplyr::select(id, gene_short_name, classification) %>%
  distinct()

joint_matrix_aggregated %>%
  filter(dose == "0",
         classification == "protein_coding") %>%
  ggplot(aes(x = log10(bulk), y = log10(sci))) +
  geom_abline(color = "black", size = .1) +
  geom_point(stroke = 0,
             size = .25,
             alpha = .25) +
  geom_smooth(method = "lm", size = .5, se = F) +
  stat_cor(
    size = 2,
    label.y.npc = 0.01,
    color = "grey31",
    method = "pearson"
  ) +
  facet_wrap( ~ cell_type + dose, ncol = 2) +
  theme(text = element_text(size = 6),
        strip.text.x = element_text(size = 6)) +
  monocle_theme_opts() +
  xlab("bulk-RNA-seq log10(TPM + 1)") +
  ylab("sci-RNA-seq log10(TPM + 1)") +
  ggsave(
    "supplemental_figure9_vehicle_correlation_per_gene_counts_prot.png",
    width = 2.25,
    height = 1.5,
    dpi = 300
  )


# Run DeSeq2

run_deseq2 = function(coldata, count_matrix) {
  coldata =
    coldata %>%
    as.data.frame()
  rownames(coldata) = coldata$sample
  
  coldata$dose = as.numeric(coldata$dose)
  #coldata$dose = relevel(coldata$dose, ref  = "0")
  
  samples = coldata$sample
  this.count_matrix = count_matrix[, samples]
  
  dds <- DESeqDataSetFromMatrix(
    countData = this.count_matrix,
    colData = coldata,
    design = ~ log(dose + .01)
  )
  dds <- DESeq(dds)
  dds
}

cell_groups =
  coldata %>%
  filter(catalog_number != "S0000") %>%
  mutate(cell_line = cell_type) %>%
  group_by(cell_type, catalog_number) %>%
  nest() %>%
  mutate(samples = purrr::map(
    data,
    .f = function(coldata_subset, coldata) {
      this_cell_type = coldata_subset$cell_line
      vehicle_samples = coldata %>%
        filter(cell_type == this_cell_type,
               catalog_number == "S0000") %>%
        pull(cell)
      sample = c(coldata_subset$cell, vehicle_samples)
      tibble(sample)
    },
    coldata
  )) %>%
  dplyr::select(samples) %>%
  mutate(sample_id = row_number()) %>%
  unnest(samples) %>%
  left_join(coldata, by = c("sample" = "cell"))

deseq_res =
  cell_groups %>%
  group_by(sample_id) %>%
  nest() %>%
  mutate(res = purrr::map(data, .f = run_deseq2, count_mat_merged))

deseq_res =
  deseq_res %>%
  mutate(res_names = purrr::map(res, resultsNames))

deseq_res_names =
  deseq_res %>%
  unnest(res_names)

deseq_res =
  deseq_res %>%
  dplyr::select(-res_names) %>%
  left_join(deseq_res_names, by = "sample_id")


get_deseq_res = function(dds_res, name_contrast) {
  DESeq2::results(dds_res, name = name_contrast)
}


deseq_res =
  deseq_res %>%
  mutate(res_df = purrr::map2(res, res_names, get_deseq_res))

deseq_to_tibble = function(deseq_res) {
  tibble(
    gene_id = deseq_res@rownames,
    base_mean = deseq_res[["baseMean"]],
    log_2_fold_change = deseq_res[["log2FoldChange"]],
    lfc_standard_error = deseq_res[["lfcSE"]],
    stat = deseq_res[["stat"]],
    p_value = deseq_res[["pvalue"]],
    p_adj = deseq_res[["padj"]]
  )
}


deseq_res =
  deseq_res %>%
  mutate(res_df = purrr::map(res_df, deseq_to_tibble))

sample_ids =
  cell_groups %>%
  filter(catalog_number != "S0000") %>%
  dplyr::select(sample_id, catalog_number, cell_type)

deseq_res =
  deseq_res %>%
  dplyr::select(sample_id, res_names, res_df) %>%
  left_join(sample_ids, by = "sample_id") %>%
  unnest() %>%
  mutate(term = paste(catalog_number, cell_type, res_names, sep = "_"))


write.table(
  deseq_res,
  "deseq_res.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)



A549_test_res = read.table("Supplementary_Table_5_A549.txt",
                           header = T)

A549_test_res =
  A549_test_res %>%
  filter(treatment %in% deseq_res$catalog_number) %>%
  filter(!grepl("Intercept", term))

A549_deseq_res =
  deseq_res %>%
  filter(!grepl("Intercept", term),
         cell_type == "A549")


A549_joint_res =
  left_join(A549_test_res,
            A549_deseq_res,
            c("id" = "gene_id",
              "treatment" = "catalog_number")) %>%
  distinct()

A549_joint_res =
  A549_joint_res %>%
  filter(p_value.x < .05 | p_value.y < .05)


A549_joint_res %>%
  drop_na() %>%
  group_by(treatment) %>%
  nest() %>%
  mutate(correlation = purrr::map(
    data,
    .f = function(data_subset) {
      cor(data_subset$estimate, data_subset$log_2_fold_change)
    }
  )) %>%
  dplyr::select(treatment, correlation) %>%
  unnest()


K562_test_res = read.table("Supplementary_Table_5_K562.txt",
                           header = T)

K562_test_res =
  K562_test_res %>%
  filter(treatment %in% deseq_res$catalog_number) %>%
  filter(!grepl("Intercept", term))


K562_deseq_res =
  deseq_res %>%
  filter(!grepl("Intercept", term),
         cell_type == "K562")


K562_joint_res =
  left_join(K562_test_res,
            K562_deseq_res,
            c("id" = "gene_id",
              "treatment" = "catalog_number")) %>%
  distinct()

K562_joint_res =
  K562_joint_res %>%
  filter(p_value.x < .01 | p_value.y < .01)

K562_joint_res %>%
  filter(p_value.x < .01 | p_value.y < .01) %>%
  group_by(treatment, cell_type.x) %>%
  summarise(n = n())


joint_res = rbind(A549_joint_res, K562_joint_res)

# Pull in drug annotations to annotate DEG results
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

joint_res = left_join(joint_res,
                      drug_annotations,
                      by = c("treatment" = "catalog_number"))
joint_res$product_name = joint_res$product_name %>% as.character()

joint_res$product_name[joint_res$product_name == "Patupilone (EPO906, Epothilone B)"] <-
  "Patupilone"

ggplot(
  joint_res %>%
    filter(p_value.x < .01, p_value.y < .01),
  aes(x = estimate, y = log_2_fold_change)
) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_vline(xintercept = 0, size = 0.25) +
  geom_point(stroke = 0, size = .25) +
  geom_abline(size = .1) +
  geom_smooth(method = "lm", size = .5, se = F) +
  stat_cor(
    size = 2,
    color = "grey31",
    label.x.npc = "left",
    label.y.npc = "top"
  ) +
  facet_grid(cell_type.x ~ product_name) +
  monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ylab("DeSeq2 Estimate") +
  xlab("sciPlex Estimate Beta") +
  ggsave(
    "supplemental_figure9_correlations_effect_sizes_.pdf",
    height = 2,
    width = 6
  )


joint_res %>%
  filter(p_value.x < 0.01, p_value.y < .01) %>%
  mutate(fold_diff = estimate / log_2_fold_change) %>%
  arrange(desc(abs(fold_diff))) %>%
  ggplot() +
  geom_density(aes(x = log2(abs(fold_diff))), fill = "grey80") +
  monocle_theme_opts() +
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6)) +
  facet_grid(cell_type.x ~ treatment) +
  xlab("log10(Absolute ratio \nof estimates)") +
  ggsave(
    "supplemental_figure9_fold_change_estimates.pdf",
    height = 2,
    width = 5
  )

joint_res %>%
  mutate(fold_diff = estimate / log_2_fold_change) %>%
  arrange(desc(abs(fold_diff))) %>%
  ggplot() +
  geom_point(
    aes(x = log10(abs(fold_diff)), y = -log10(p_value.x)),
    color = "grey80",
    size = .25,
    alpha = .5
  ) +
  monocle_theme_opts() +
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6)) +
  geom_vline(xintercept = 1,
             color = "red",
             size = 0.5) +
  geom_vline(xintercept = -1,
             color = "red",
             size = 0.5) +
  
  xlab("log10(Absolute ratio \nof estimates)") +
  ggsave(
    "supplemental_figure9_fold_change_estimates_scatter.pdf",
    height = 1.5,
    width = 1.5
  )
