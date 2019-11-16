# Need to download data from  GSE92742 for L1000 comparison
# Data not included in repository


# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level1_LXB_n1403502.tar.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level2_GEX_delta_n49216x978.gctx.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level2_GEX_epsilon_n1269922x978.gctx.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level4_ZSPCINF_mlr12k_n1319138x12328.gctx.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_README.pdf
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_auxiliary_datasets.tar.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_inst_info.txt.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_metrics.txt.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_metrics.txt.gz
# !Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_SHA512SUMS.txt.gz


# devtools::install_github(repo = 'cole-trapnell-lab/monocle3', ref = "d4a9a35")

# Input path to github directory
path_to_github = "~/sci-Plex/"


# Set directory for sciPlex github 
bin_directory = paste0(path_to_github,
                       "bin",
                       sep = "")

suppressPackageStartupMessages({
  library(cmapR)
  library(stringr)
  library(dplyr)
  library(pheatmap)
  library(ggplot2)
  library(tidyr)
  library(ggrepel)
  library(Matrix)
  library(broom)
  library(purrr)
  library(ggpubr)
  library(monocle3)
})

# Set DelayedArray Parameters
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size = 1000e7)


compute_effect_matrix = function(test_res) {
  dose_terms =
    test_res %>%
    dplyr::filter(grepl("dose_", term))
  
  effect_matrix =
    dose_terms %>%
    ungroup() %>%
    dplyr::select(id, treatment, estimate) %>%
    spread(treatment, estimate, fill = 0)
  
  colnames(effect_matrix) = c("id",  colnames(effect_matrix)[-1])
  
  gene_ids = as.character(effect_matrix$id)
  effect_matrix = effect_matrix[,-1]
  effect_matrix = as.matrix(effect_matrix)
  effect_matrix[is.na(effect_matrix)] = 0
  row.names(effect_matrix) = gene_ids
  
  return(effect_matrix)
}

# Read in sciPlex DEGs

A549_test_res = read.table("Supplementary_Table_5_A549.txt",
                           header = T)
MCF7_test_res = read.table("Supplementary_Table_5_MCF7.txt",
                           header = T)

A549_effect_matrix = compute_effect_matrix(A549_test_res)

MCF7_effect_matrix = compute_effect_matrix(MCF7_test_res)

gene_rowdata =
  data.frame(id = A549_test_res$id,
             gene_short_name = A549_test_res$gene_short_name) %>%
  distinct()

gene_rowdata =
  gene_rowdata %>%
  filter(id %in% row.names(A549_effect_matrix))


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

######################


# Level 5 connectivity map data
# Contains data for all genes including measured and inferred

# Download from GSE92742
ds_path <-
  "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
my_ds <- parse.gctx(ds_path)

##### Adjust column metadata #####
# Download from GSE92742
col_meta_path <-
  "GSE92742_Broad_LINCS_sig_info.txt"
col_meta <-
  read.delim(col_meta_path, sep = "\t", stringsAsFactors = F)

# get column metadata,
col_meta_from_gctx <- read.gctx.meta(ds_path, dim = "col")

col_meta =
  col_meta %>%
  dplyr::mutate(sample_count = stringr::str_count(string = distil_id, pattern = "\\|") + 1)

max_samples = max(col_meta$sample_count)

distil_id_samples = as.character(c())
for (i in seq(1, max_samples)) {
  new_sample_col_id = as.character(paste("id", i, sep = ""))
  distil_id_samples = c(distil_id_samples, new_sample_col_id)
}

col_meta =
  col_meta %>%
  dplyr::select(-sample_count)

sample_keys =
  col_meta %>%
  tidyr::separate(distil_id, into = distil_id_samples, sep = "\\|") %>%
  dplyr::select(sig_id, dplyr::starts_with("id"))

sample_keys =
  sample_keys %>%
  tidyr::gather(starts_with("id"),
                key = "dummy",
                value = "id",
                na.rm = T) %>%
  dplyr::select(-dummy) %>%
  distinct() %>%
  dplyr::filter(!grepl(pattern = ":-666", x = id))

col_meta =
  col_meta %>%
  dplyr::select(-distil_id)

col_meta_long = dplyr::left_join(sample_keys, col_meta, by = "sig_id")

# remove duplicated rows corresponding to different sig_ids
# drop rows with missing metadata
col_meta_long =
  col_meta_long %>%
  dplyr::distinct()

# there are still duplicated id rows which correspond to duplicated perturbation ids
col_meta_long_duplicated =
  col_meta_long %>%
  dplyr::group_by(sig_id) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::ungroup()

col_meta_long_duplicated =
  col_meta_long_duplicated %>%
  dplyr::group_by(sig_id) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup()

col_meta_long =
  col_meta_long %>%
  dplyr::group_by(sig_id) %>%
  dplyr::filter(dplyr::n() == 1) %>%
  dplyr::ungroup()

col_meta_long = rbind(col_meta_long, col_meta_long_duplicated)

num_ids =
  col_meta_long$sig_id %>%
  length()

unique_ids =
  col_meta_long$sig_id %>%
  unique() %>%
  length()

num_ids == unique_ids

col_meta_from_gctx %>%
  dplyr::mutate(meta_id_present = (id %in% col_meta_long$sig_id)) %>%
  group_by(meta_id_present) %>%
  summarise(n = n())

col_meta_from_gctx %>%
  mutate(meta_id_present = (id %in% col_meta_long$sig_id)) %>%
  dplyr::filter(!meta_id_present) %>%
  tail()

col_meta_from_gctx =
  left_join(col_meta_from_gctx, col_meta_long, by = c("id" = "sig_id"))


selleckchem_to_l1000 = read.table(
  file = "bin/sciPlex_l1000_mapping.txt",
  sep = "\t",
  col.names = c("catalog_number", "pert_id")
)

# Add the selleckPlex catalog numbers to the l1000 meta_data
col_meta_from_gctx =
  left_join(col_meta_from_gctx, selleckchem_to_l1000, by = "pert_id")

# Download from GSE92742
# Load Level 2 data

ds_path_level_2 <-
  "GSE92742_Broad_LINCS_Level2_GEX_epsilon_n1269922x978.gctx"
my_ds_level_2 <- parse.gctx(ds_path_level_2)
genes_measured_1000 = rownames(my_ds_level_2@mat)


####### Adjust row metadata #######

# Download from GSE92742
row_meta_path <-
  "GSE92742_Broad_LINCS_gene_info.txt"
row_meta <-
  read.delim(row_meta_path, sep = "\t", stringsAsFactors = F)

row_meta =
  row_meta %>%
  dplyr::select(everything(), gene_short_name = pr_gene_symbol) %>%
  dplyr::mutate(pr_gene_id = as.character(pr_gene_id))

row_meta =
  left_join(row_meta, gene_rowdata, by = "gene_short_name") %>%
  dplyr::filter(!is.na(id))

# Remove genes that dont have a matching ensmbl gene id

row_meta =
  row_meta %>%
  dplyr::group_by(pr_gene_id) %>%
  dplyr::filter(n() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(id = as.character(id)) %>%
  dplyr::mutate(measured_l1000 = pr_gene_id %in% genes_measured_1000)

row_meta = row_meta[row_meta$pr_gene_id %in% rownames(my_ds@mat),]

####### Get A549 data from Connectivity Map data #######

a549_data_columns =
  col_meta_from_gctx %>%
  dplyr::filter(
    pert_id %in% selleckchem_to_l1000$pert_id,
    cell_id == "A549",
    pert_time == 24,
    pert_time_unit == "h"
  ) %>%
  dplyr::pull(id)

a549_data_mat = my_ds@mat[, (colnames(my_ds@mat) %in% a549_data_columns)]

a549_data_mat = a549_data_mat[row_meta$pr_gene_id,]

identical(row.names(a549_data_mat), row_meta$pr_gene_id)

A549_effect_matrix = A549_effect_matrix[rownames(A549_effect_matrix) %in% row_meta$id,]
identical(row.names(A549_effect_matrix), row_meta$id)


####### Get MCF7 data from Connectivity Map data #######
mcf7_data_columns =
  col_meta_from_gctx %>%
  dplyr::filter(
    pert_id %in% selleckchem_to_l1000$pert_id,
    cell_id == "MCF7",
    pert_time == 24,
    pert_time_unit == "h"
  ) %>%
  dplyr::pull(id)
mcf7_data_mat = my_ds@mat[, (colnames(my_ds@mat) %in% mcf7_data_columns)]

mcf7_data_mat = mcf7_data_mat[row_meta$pr_gene_id,]
identical(row.names(mcf7_data_mat), row_meta$pr_gene_id)


MCF7_effect_matrix = MCF7_effect_matrix[rownames(MCF7_effect_matrix) %in% row_meta$id,]
identical(row.names(A549_effect_matrix), row_meta$id)

####### Make plots for comparisons #######

l1000_sciPlex_corr = function(coldata_subset,
                              row_meta,
                              l1000_coldata,
                              effect_matrix,
                              test_res,
                              l1000_data_mat) {
  querry_catalog_number =
    coldata_subset$catalog_number %>%
    unique() %>%
    as.character()
  
  print(querry_catalog_number)
  
  sciPlex_df = data.frame(betas_sciPlex = effect_matrix[, querry_catalog_number] %>% as.numeric,
                          id = effect_matrix[, querry_catalog_number] %>% names)
  
  sciPlex_DEGs =
    test_res %>%
    filter(grepl(querry_catalog_number, term)) %>%
    dplyr::select(id, term, q_value, p_value)
  
  sciPlex_df =
    sciPlex_df %>%
    left_join(sciPlex_DEGs, by = "id") %>%
    left_join(row_meta, by = "id") %>%
    dplyr::select(betas_sciPlex,
                  ensmble_id = id,
                  term,
                  q_value,
                  p_value,
                  pr_gene_id) %>%
    drop_na()
  
  colnames(l1000_data_mat)
  
  if (sum(colnames(l1000_data_mat) %in% coldata_subset$id) == 1) {
    tmp.vector  =  l1000_data_mat[, colnames(l1000_data_mat) %in% coldata_subset$id]
    l1000_df = data.frame(
      id = colnames(l1000_data_mat)[colnames(l1000_data_mat) %in% coldata_subset$id],
      l1000_modz = tmp.vector %>% as.numeric,
      pr_gene_id = names(tmp.vector)
    )
  } else{
    l1000_data_mat = l1000_data_mat[, colnames(l1000_data_mat) %in% coldata_subset$id]
    
    l1000_df =
      l1000_data_mat %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "pr_gene_id") %>%
      gather(key = id, value = l1000_modz, -pr_gene_id)
  }
  
  l1000_df =
    l1000_df %>%
    left_join(coldata_subset, by = "id") %>%
    dplyr::select(everything(), -id.y) %>%
    left_join(row_meta %>% dplyr::select(everything(), ensmble_id = id), by = "pr_gene_id")
  
  
  l1000_df =
    l1000_df %>%
    group_by(pert_idose, gene_short_name) %>%
    mutate(mean_l1000_signal = mean(l1000_modz),
           std_l1000 = sqrt(var(l1000_modz))) %>%
    ungroup() %>%
    mutate(pert_idose_numeric = stringr::str_split_fixed(pert_idose, " ", 2)[, 1])
  
  joint_df =
    left_join(l1000_df, sciPlex_df)
  
  max_concentration =
    joint_df %>%
    filter(pert_dose_unit == "µM") %>%
    pull(pert_idose_numeric) %>%
    as.numeric %>%
    max()
  
  df_max =
    joint_df %>%
    filter(
      as.numeric(pert_idose_numeric) == max_concentration,
      pert_dose_unit == "µM",
      measured_l1000,
      p_value < .01
    )
  
  return(
    cor.test(df_max$mean_l1000_signal, df_max$betas_sciPlex, method = "spearman") %>%
      tidy()
  )
}

correlations_mcf7 =
  col_meta_from_gctx %>%
  filter(
    pert_id %in% selleckchem_to_l1000$pert_id,
    cell_id == "MCF7",
    pert_time == 24,
    pert_time_unit == "h",
    !is.na(catalog_number)
  ) %>%
  mutate(catalog_dummy = catalog_number) %>%
  group_by(catalog_dummy) %>%
  nest() %>%
  mutate(
    correlation = purrr::map(
      data,
      .f = l1000_sciPlex_corr,
      row_meta,
      col_meta_from_gctx,
      MCF7_effect_matrix,
      MCF7_test_res,
      mcf7_data_mat
    )
  )  %>%
  dplyr::select(-data) %>%
  unnest() %>%
  dplyr::rename(catalog_number = catalog_dummy) %>%
  left_join(drug_annotations, by = "catalog_number")



correlations_a549 =
  col_meta_from_gctx %>%
  dplyr::filter(
    pert_id %in% selleckchem_to_l1000$pert_id,
    cell_id == "A549",
    pert_time == 24,
    pert_time_unit == "h",
    !is.na(catalog_number)
  ) %>%
  dplyr::mutate(catalog_dummy = catalog_number) %>%
  dplyr::group_by(catalog_dummy) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    correlation = purrr::map(
      data,
      .f = l1000_sciPlex_corr,
      row_meta,
      col_meta_from_gctx,
      A549_effect_matrix,
      A549_test_res,
      a549_data_mat
    )
  )  %>%
  dplyr::select(-data) %>%
  unnest() %>%
  dplyr::rename(catalog_number = catalog_dummy) %>%
  left_join(drug_annotations, by = "catalog_number")


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


correlations_mcf7 %>%
  filter(p.value < .05) %>%
  group_by(pathway_level_1) %>%
  mutate(median = median(estimate)) %>%
  ggplot(aes(x = reorder(pathway_level_1, median), y = estimate)) +
  geom_boxplot(
    aes(fill = pathway_level_1),
    color = "black",
    size = .25,
    outlier.size = 0,
    outlier.stroke = 0
  ) +
  # geom_jitter(data = correlations_a549 %>% filter(catalog_number != "S1191"),
  #             color = "black", size = 1, stroke = 0) +
  # geom_jitter(data = correlations_a549 %>% filter(catalog_number == "S1191"),
  #             color = "red", size = 1, stroke = 0) +
  scale_fill_manual("Pathway", values = pathway_level_1_colors) +  coord_flip() +
  monocle_theme_opts() +
  theme(legend.position = "none",
        text = element_text(size = 6)) +
  xlab("") +
  ylab("Spearman Correlation") +
  ggsave(
    "supplemental_figure10_correlations_all_mcf7.pdf",
    height = 2.5,
    width = 2.5
  )


correlations_mcf7_medians =
  correlations_mcf7 %>%
  filter(p.value < .05) %>%
  group_by(pathway_level_1) %>%
  mutate(median = median(estimate))

ggplot(correlations_mcf7_medians, aes(x = reorder(pathway_level_1, median), y = estimate)) +
  geom_boxplot(
    aes(fill = pathway_level_1),
    color = "black",
    size = .25,
    outlier.size = 0,
    outlier.stroke = 0
  ) +
  geom_jitter(aes(color = (catalog_number != "S1191")), size = 1, stroke = 0) +
  scale_fill_manual("Pathway", values = pathway_level_1_colors) +
  scale_color_manual(values = c("red", "black")) +
  coord_flip() +
  monocle_theme_opts() +
  theme(legend.position = "none",
        text = element_text(size = 6)) +
  xlab("") +
  ylim(-.5, 0.75) +
  ylab("Spearman Correlation") +
  ggsave(
    "supplemental_figure10_correlations_all_mcf7.pdf",
    height = 2.5,
    width = 2.5
  )

correlations_a549_medians =
  correlations_a549 %>%
  filter(p.value < .05) %>%
  group_by(pathway_level_1) %>%
  mutate(median = median(estimate))

ggplot(correlations_a549_medians, aes(x = reorder(pathway_level_1, median), y = estimate)) +
  geom_boxplot(
    aes(fill = pathway_level_1),
    color = "black",
    size = .25,
    outlier.size = 0,
    outlier.stroke = 0
  ) +
  geom_jitter(aes(color = (catalog_number != "S1628")), size = 1, stroke = 0) +
  scale_fill_manual("Pathway", values = pathway_level_1_colors) +
  scale_color_manual(values = c("red", "black")) +
  coord_flip() +
  monocle_theme_opts() +
  theme(legend.position = "none",
        text = element_text(size = 6)) +
  xlab("") +
  ylab("Spearman Correlation") +
  ggsave(
    "supplemental_figure10_correlations_all_a549.pdf",
    height = 2.5,
    width = 2.5
  )


###############


l1000_sciPlex_match = function(coldata_subset,
                               row_meta,
                               l1000_coldata,
                               effect_matrix,
                               test_res,
                               l1000_data_mat) {
  querry_catalog_number =
    coldata_subset$catalog_number %>%
    unique() %>%
    as.character()
  
  print(querry_catalog_number)
  
  sciPlex_df = data.frame(betas_sciPlex = effect_matrix[, querry_catalog_number] %>% as.numeric,
                          id = effect_matrix[, querry_catalog_number] %>% names)
  
  sciPlex_DEGs =
    test_res %>%
    filter(grepl(querry_catalog_number, term)) %>%
    dplyr::select(id, term, q_value, p_value)
  
  sciPlex_df =
    sciPlex_df %>%
    left_join(sciPlex_DEGs, by = "id") %>%
    left_join(row_meta, by = "id") %>%
    dplyr::select(betas_sciPlex,
                  ensmble_id = id,
                  term,
                  q_value,
                  p_value,
                  pr_gene_id) %>%
    drop_na()
  
  colnames(l1000_data_mat)
  
  if (sum(colnames(l1000_data_mat) %in% coldata_subset$id) == 1) {
    tmp.vector  =  l1000_data_mat[, colnames(l1000_data_mat) %in% coldata_subset$id]
    l1000_df = data.frame(
      id = colnames(l1000_data_mat)[colnames(l1000_data_mat) %in% coldata_subset$id],
      l1000_modz = tmp.vector %>% as.numeric,
      pr_gene_id = names(tmp.vector)
    )
  } else{
    l1000_data_mat = l1000_data_mat[, colnames(l1000_data_mat) %in% coldata_subset$id]
    
    l1000_df =
      l1000_data_mat %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "pr_gene_id") %>%
      gather(key = id, value = l1000_modz, -pr_gene_id)
  }
  
  l1000_df =
    l1000_df %>%
    left_join(coldata_subset, by = "id") %>%
    dplyr::select(everything(), -id.y) %>%
    left_join(row_meta %>% dplyr::select(everything(), ensmble_id = id), by = "pr_gene_id")
  
  
  l1000_df =
    l1000_df %>%
    group_by(pert_idose, gene_short_name) %>%
    mutate(mean_l1000_signal = mean(l1000_modz),
           std_l1000 = sqrt(var(l1000_modz))) %>%
    ungroup() %>%
    mutate(pert_idose_numeric = stringr::str_split_fixed(pert_idose, " ", 2)[, 1])
  
  joint_df =
    left_join(l1000_df, sciPlex_df)
  
  return(joint_df)
}





mcf7_l1000_sciPlex_match =
  col_meta_from_gctx %>%
  filter(
    pert_id %in% selleckchem_to_l1000$pert_id,
    cell_id == "MCF7",
    pert_time == 24,
    pert_time_unit == "h",
    !is.na(catalog_number)
  ) %>%
  mutate(catalog_dummy = catalog_number) %>%
  filter(catalog_number %in% c("S1045", "S1191", "S1628")) %>%
  group_by(catalog_dummy) %>%
  nest() %>%
  mutate(
    correlation = purrr::map(
      data,
      .f = l1000_sciPlex_match,
      row_meta,
      col_meta_from_gctx,
      MCF7_effect_matrix,
      MCF7_test_res,
      mcf7_data_mat
    )
  )  %>%
  dplyr::select(-data) %>%
  unnest()


mcf7_l1000_sciPlex_match %>%
  filter(q_value < .05,
         catalog_dummy == "S1045",
         measured_l1000) %>%
  filter(pert_idose %in% c("10 µM", "1 µM", "500 nM", "100 nM",
                           "5 µM" , "10 nM", "1 nM", "0.04 µM")) %>%
  mutate(pert_idose_factor =  factor(
    pert_idose,
    levels = c("1 nM", "10 nM", "0.04 µM", "100 nM", "500 nM", "1 µM", "5 µM", "10 µM")
  )) %>%
  dplyr::select(betas_sciPlex,
                mean_l1000_signal,
                pert_idose,
                pert_idose_factor) %>%
  distinct() %>%
  ggplot(aes(x = betas_sciPlex, y = mean_l1000_signal)) +
  geom_smooth(method = lm, size = .5, se = F) +
  geom_point(color = "black",
             size = .85,
             stroke = 0) +
  geom_point(aes(color = pert_idose_factor),
             size = .55,
             stroke = 0) +
  stat_cor(size = 2,
           label.y.npc = 0.01,
           color = "grey31") +
  monocle3::monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    strip.text = element_text(size = 6)
  ) +
  facet_wrap( ~ pert_idose_factor, ncol = 2) +
  scale_color_viridis_d(option = "magma") +
  xlab("sciPlex Dose Dependent Beta") +
  ylab("L1000 MODZ Score") +
  ggsave(
    "supplemental_figure10_l1000_measured_tsa_mcf7.pdf",
    height = 5,
    width = 2.5
  )


mcf7_l1000_sciPlex_match %>%
  filter(
    q_value < .05,
    catalog_dummy == "S1191",
    measured_l1000,
    pert_idose_numeric == 10,
    pert_dose_unit == "µM"
  )  %>%
  dplyr::select(betas_sciPlex,
                mean_l1000_signal,
                pert_idose) %>%
  distinct() %>%
  ggplot(aes(x = betas_sciPlex, y = mean_l1000_signal)) +
  geom_smooth(method = lm, size = .5, se = F) +
  geom_point(color = "black",
             size = .55,
             stroke = 0) +
  stat_cor(size = 2,
           label.y.npc = 0.01,
           color = "grey31") +
  geom_point(color = "chartreuse3",
             size = .45,
             stroke = 0) +
  monocle3::monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    strip.text = element_text(size = 6)
  ) +
  xlab("sciPlex Dose Beta") +
  ylab("L1000 ModZ Score") +
  ggsave(
    "supplemental_figure10_l1000_measured_fuv_10_mcf7.pdf",
    height = 1.5,
    width = 1.5
  )


mcf7_l1000_sciPlex_match %>%
  filter(
    q_value < .05,
    catalog_dummy == "S1628",
    measured_l1000,
    pert_idose_numeric == 10,
    pert_dose_unit == "µM"
  )  %>%
  dplyr::select(betas_sciPlex,
                mean_l1000_signal,
                pert_idose) %>%
  distinct() %>%
  ggplot(aes(x = betas_sciPlex, y = mean_l1000_signal)) +
  geom_smooth(method = lm, size = .5, se = F) +
  geom_point(color = "black",
             size = .55,
             stroke = 0) +
  stat_cor(size = 2,
           label.y.npc = 0.01,
           color = "grey31") +
  geom_point(color = "chartreuse3",
             size = .45,
             stroke = 0) +
  monocle3::monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    strip.text = element_text(size = 6)
  ) +
  xlab("sciPlex Dose Beta") +
  ylab("L1000 ModZ Score") +
  ggsave(
    "supplemental_figure10_l1000_measured_triam_10_mcf7.pdf",
    height = 1.5,
    width = 1.5
  )




#######  A549 #######


a549_l1000_sciPlex_match =
  col_meta_from_gctx %>%
  filter(
    pert_id %in% selleckchem_to_l1000$pert_id,
    cell_id == "A549",
    pert_time == 24,
    pert_time_unit == "h",
    !is.na(catalog_number)
  ) %>%
  mutate(catalog_dummy = catalog_number) %>%
  filter(catalog_number %in% c("S1045", "S1191", "S1628")) %>%
  group_by(catalog_dummy) %>%
  nest() %>%
  mutate(
    correlation = purrr::map(
      data,
      .f = l1000_sciPlex_match,
      row_meta,
      col_meta_from_gctx,
      A549_effect_matrix,
      A549_test_res,
      a549_data_mat
    )
  )  %>%
  dplyr::select(-data) %>%
  unnest()



a549_l1000_sciPlex_match %>%
  filter(
    q_value < .05,
    catalog_dummy == "S1191",
    measured_l1000,
    pert_idose_numeric == 10,
    pert_dose_unit == "µM"
  )  %>%
  dplyr::select(betas_sciPlex,
                mean_l1000_signal,
                pert_idose) %>%
  distinct() %>%
  ggplot(aes(x = betas_sciPlex, y = mean_l1000_signal)) +
  geom_smooth(method = lm, size = .5, se = F) +
  geom_point(color = "black",
             size = .55,
             stroke = 0) +
  stat_cor(size = 2,
           label.y.npc = 0.01,
           color = "grey31") +
  geom_point(color = "chartreuse3",
             size = .45,
             stroke = 0) +
  monocle3::monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    strip.text = element_text(size = 6)
  ) +
  xlab("sciPlex Dose Beta") +
  ylab("L1000 ModZ Score") +
  ggsave(
    "supplemental_figure10_l1000_measured_fuv_10_a549.pdf",
    height = 1.5,
    width = 1.5
  )


a549_l1000_sciPlex_match %>%
  filter(
    q_value < .05,
    catalog_dummy == "S1628",
    measured_l1000,
    pert_idose_numeric == 10,
    pert_dose_unit == "µM"
  )  %>%
  dplyr::select(betas_sciPlex,
                mean_l1000_signal,
                pert_idose) %>%
  distinct() %>%
  ggplot(aes(x = betas_sciPlex, y = mean_l1000_signal)) +
  geom_smooth(method = lm, size = .5, se = F) +
  geom_point(color = "black",
             size = .55,
             stroke = 0) +
  stat_cor(size = 2,
           label.y.npc = 0.01,
           color = "grey31") +
  geom_point(color = "chartreuse3",
             size = .45,
             stroke = 0) +
  monocle3::monocle_theme_opts() +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    strip.text = element_text(size = 6)
  ) +
  xlab("sciPlex Dose Beta") +
  ylab("L1000 ModZ Score") +
  ggsave(
    "supplemental_figure10_l1000_measured_triam_10_a549.pdf",
    height = 1.5,
    width = 1.5
  )
