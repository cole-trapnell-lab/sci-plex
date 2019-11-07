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
  library(tidymodels)
  library(furrr)
  library(tictoc)
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
  source(paste0(bin_directory,
                "test_drug_dependence.R",
                sep = ""))
})



DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size = 1000e7)


# Set DelayedArray Parameters
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size = 1000e7)

# User defined functions

setwd(paste0(path_to_github,
             "large_screen",
             sep = ""))

path_to_cds = "MCF7_24hrs.RDS"

cds = readRDS(path_to_cds)
cds = estimate_cell_cycle(cds, cc.genes$s.gene, cc.genes$g2m.gene)


vehicle_cells = colData(cds) %>% as.data.frame() %>% filter(vehicle == TRUE)
vehicle_cells = as.character(vehicle_cells$cell)


options(future.globals.maxSize = 891289600 * 10)
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)
plan(multiprocess)
tic()
adjusted_model_test_res =
  colData(cds) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  mutate(
    test_res = furrr::future_map(
      data,
      .f = test_drug_dependence,
      cds,
      vehicle_cells,
      residualModelFormulaTerms = NULL,
      .progress = TRUE
    )
  )
adjusted_model_test_res = adjusted_model_test_res %>% dplyr::select(-data) %>% unnest()
toc()


# Correct for multiple testing
adjusted_model_test_res = adjusted_model_test_res %>% group_by(term) %>% mutate(q_value = p.adjust(p_value))

write.table(
  adjusted_model_test_res,
  "Supplementary_Table_5_MCF7.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)




low_PI_cells =
  cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(proliferation_index <= 2.25) %>%
  pull(cell)

cds_low_proliferation = cds[, low_PI_cells]



vehicle_cells =
  colData(cds_low_proliferation) %>%
  as.data.frame() %>%
  filter(vehicle == TRUE) %>%
  pull(cell) %>%
  as.character()



options(future.globals.maxSize = 891289600 * 10)
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)
plan(multiprocess)
tic()
adjusted_model_test_res =
  colData(cds_low_proliferation) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  mutate(
    test_res = furrr::future_map(
      data,
      .f = test_drug_dependence,
      cds_low_proliferation,
      vehicle_cells,
      residualModelFormulaTerms = NULL,
      .progress = TRUE
    )
  )
adjusted_model_test_res = adjusted_model_test_res %>% dplyr::select(-data) %>% unnest()
toc()


# Correct for multiple testing
adjusted_model_test_res = adjusted_model_test_res %>% group_by(term) %>% mutate(q_value = p.adjust(p_value))

write.table(
  adjusted_model_test_res,
  "MCF7_low_proliferation_adjusted_model_test_res.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)



high_PI_cells =
  cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(proliferation_index > 2.25) %>%
  pull(cell)

cds_high_proliferation = cds[, high_PI_cells]



vehicle_cells =
  colData(cds_high_proliferation) %>%
  as.data.frame() %>%
  filter(vehicle == TRUE) %>%
  pull(cell) %>%
  as.character()



options(future.globals.maxSize = 891289600 * 10)
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)
plan(multiprocess)
tic()
adjusted_model_test_res =
  colData(cds_high_proliferation) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  mutate(
    test_res = furrr::future_map(
      data,
      .f = test_drug_dependence,
      cds_high_proliferation,
      vehicle_cells,
      residualModelFormulaTerms = NULL,
      .progress = TRUE
    )
  )
adjusted_model_test_res = adjusted_model_test_res %>% dplyr::select(-data) %>% unnest()
toc()


# Correct for multiple testing
adjusted_model_test_res = adjusted_model_test_res %>% group_by(term) %>% mutate(q_value = p.adjust(p_value))

write.table(
  adjusted_model_test_res,
  "MCF7_high_proliferation_adjusted_model_test_res.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)

#################### A549 #################### 


path_to_cds = "A549_24hrs.RDS"

cds = readRDS(path_to_cds)
cds = estimate_cell_cycle(cds, cc.genes$s.gene, cc.genes$g2m.gene)


vehicle_cells = colData(cds) %>% as.data.frame() %>% filter(vehicle == TRUE)
vehicle_cells = as.character(vehicle_cells$cell)


options(future.globals.maxSize = 891289600 * 10)
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)
plan(multiprocess)
tic()
adjusted_model_test_res =
  colData(cds) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  mutate(
    test_res = furrr::future_map(
      data,
      .f = test_drug_dependence,
      cds,
      vehicle_cells,
      residualModelFormulaTerms = NULL,
      .progress = TRUE
    )
  )
adjusted_model_test_res = adjusted_model_test_res %>% dplyr::select(-data) %>% unnest()
toc()


# Correct for multiple testing
adjusted_model_test_res = adjusted_model_test_res %>% group_by(term) %>% mutate(q_value = p.adjust(p_value))

write.table(
  adjusted_model_test_res,
  "Supplementary_Table_5_A549.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)




low_PI_cells =
  cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(proliferation_index <= 2) %>%
  pull(cell)

cds_low_proliferation = cds[, low_PI_cells]



vehicle_cells =
  colData(cds_low_proliferation) %>%
  as.data.frame() %>%
  filter(vehicle == TRUE) %>%
  pull(cell) %>%
  as.character()



options(future.globals.maxSize = 891289600 * 10)
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)
plan(multiprocess)
tic()
adjusted_model_test_res =
  colData(cds_low_proliferation) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  mutate(
    test_res = furrr::future_map(
      data,
      .f = test_drug_dependence,
      cds_low_proliferation,
      vehicle_cells,
      residualModelFormulaTerms = NULL,
      .progress = TRUE
    )
  )
adjusted_model_test_res = adjusted_model_test_res %>% dplyr::select(-data) %>% unnest()
toc()


# Correct for multiple testing
adjusted_model_test_res = adjusted_model_test_res %>% group_by(term) %>% mutate(q_value = p.adjust(p_value))

write.table(
  adjusted_model_test_res,
  "A549_low_proliferation_adjusted_model_test_res.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)



high_PI_cells =
  cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(proliferation_index > 2) %>%
  pull(cell)

cds_high_proliferation = cds[, high_PI_cells]



vehicle_cells =
  colData(cds_high_proliferation) %>%
  as.data.frame() %>%
  filter(vehicle == TRUE) %>%
  pull(cell) %>%
  as.character()



options(future.globals.maxSize = 891289600 * 10)
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)
plan(multiprocess)
tic()
adjusted_model_test_res =
  colData(cds_high_proliferation) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  mutate(
    test_res = furrr::future_map(
      data,
      .f = test_drug_dependence,
      cds_high_proliferation,
      vehicle_cells,
      residualModelFormulaTerms = NULL,
      .progress = TRUE
    )
  )
adjusted_model_test_res = adjusted_model_test_res %>% dplyr::select(-data) %>% unnest()
toc()


# Correct for multiple testing
adjusted_model_test_res = adjusted_model_test_res %>% group_by(term) %>% mutate(q_value = p.adjust(p_value))

write.table(
  adjusted_model_test_res,
  "A549_high_proliferation_adjusted_model_test_res.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)



####################### K562 ####################### 



path_to_cds = "K562_24hrs.RDS"

cds = readRDS(path_to_cds)
cds = estimate_cell_cycle(cds, cc.genes$s.gene, cc.genes$g2m.gene)


vehicle_cells = colData(cds) %>% as.data.frame() %>% filter(vehicle == TRUE)
vehicle_cells = as.character(vehicle_cells$cell)


options(future.globals.maxSize = 891289600 * 10)
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)
plan(multiprocess)
tic()
adjusted_model_test_res =
  colData(cds) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  mutate(
    test_res = furrr::future_map(
      data,
      .f = test_drug_dependence,
      cds,
      vehicle_cells,
      residualModelFormulaTerms = NULL,
      .progress = TRUE
    )
  )
adjusted_model_test_res = adjusted_model_test_res %>% dplyr::select(-data) %>% unnest()
toc()


# Correct for multiple testing
adjusted_model_test_res = adjusted_model_test_res %>% group_by(term) %>% mutate(q_value = p.adjust(p_value))

write.table(
  adjusted_model_test_res,
  "Supplementary_Table_5_K562.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)




low_PI_cells =
  cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(proliferation_index <= 2) %>%
  pull(cell)

cds_low_proliferation = cds[, low_PI_cells]



vehicle_cells =
  colData(cds_low_proliferation) %>%
  as.data.frame() %>%
  filter(vehicle == TRUE) %>%
  pull(cell) %>%
  as.character()



options(future.globals.maxSize = 891289600 * 10)
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)
plan(multiprocess)
tic()
adjusted_model_test_res =
  colData(cds_low_proliferation) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  mutate(
    test_res = furrr::future_map(
      data,
      .f = test_drug_dependence,
      cds_low_proliferation,
      vehicle_cells,
      residualModelFormulaTerms = NULL,
      .progress = TRUE
    )
  )
adjusted_model_test_res = adjusted_model_test_res %>% dplyr::select(-data) %>% unnest()
toc()


# Correct for multiple testing
adjusted_model_test_res = adjusted_model_test_res %>% group_by(term) %>% mutate(q_value = p.adjust(p_value))

write.table(
  adjusted_model_test_res,
  "K562_low_proliferation_adjusted_model_test_res.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)



high_PI_cells =
  cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(proliferation_index > 2) %>%
  pull(cell)

cds_high_proliferation = cds[, high_PI_cells]



vehicle_cells =
  colData(cds_high_proliferation) %>%
  as.data.frame() %>%
  filter(vehicle == TRUE) %>%
  pull(cell) %>%
  as.character()



options(future.globals.maxSize = 891289600 * 10)
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)
plan(multiprocess)
tic()
adjusted_model_test_res =
  colData(cds_high_proliferation) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  mutate(
    test_res = furrr::future_map(
      data,
      .f = test_drug_dependence,
      cds_high_proliferation,
      vehicle_cells,
      residualModelFormulaTerms = NULL,
      .progress = TRUE
    )
  )
adjusted_model_test_res = adjusted_model_test_res %>% dplyr::select(-data) %>% unnest()
toc()


# Correct for multiple testing
adjusted_model_test_res = adjusted_model_test_res %>% group_by(term) %>% mutate(q_value = p.adjust(p_value))

write.table(
  adjusted_model_test_res,
  "K562_high_proliferation_adjusted_model_test_res.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)

################### A549 72 hours ################### 

# Path to cds containing all the data
path_to_cds = "cds.RDS"

# Read in CDS and extract cells passing filters and human transcripts
cds = readRDS(path_to_cds)

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


# Isolate cells at 72 hours for DEG analysis
cells_72hrs =
  colData(cds) %>%
  as.data.frame() %>%
  filter(time_point == 72) %>%
  pull(cell)

# Subset 24 hour cells
cds_72hrs = cds[, cells_72hrs]


cds_72hrs = estimate_cell_cycle(cds_72hrs, cc.genes$s.gene, cc.genes$g2m.gene)


vehicle_cells =
  colData(cds_72hrs) %>%
  as.data.frame() %>%
  filter(vehicle == TRUE) %>%
  pull(cell) %>%
  as.character()


options(future.globals.maxSize = 891289600 * 10)
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)
plan(multiprocess)
tic()
adjusted_model_test_res =
  colData(cds_72hrs) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  mutate(
    test_res = furrr::future_map(
      data,
      .f = test_drug_dependence,
      cds_72hrs,
      vehicle_cells,
      residualModelFormulaTerms = NULL,
      .progress = TRUE
    )
  )
adjusted_model_test_res = adjusted_model_test_res %>% dplyr::select(-data) %>% unnest()
toc()


# Correct for multiple testing
adjusted_model_test_res = adjusted_model_test_res %>% group_by(term) %>% mutate(q_value = p.adjust(p_value))

write.table(
  adjusted_model_test_res,
  "Supplementary_Table_5_A549_72hours.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)


######### Testing the DEGs of vehicle proliferation high vs proliferation low cells ##############

# Path to cds containing all the data
path_to_cds = "cds.RDS"

# Read in CDS and extract cells passing filters and human transcripts
cds = readRDS(path_to_cds)

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

# Isolate cells at 24 hours for analysis
cells_24hrs =
  colData(cds) %>%
  as.data.frame() %>%
  filter(time_point == 24) %>%
  pull(cell)

cds = cds[, cells_24hrs]

cds <- estimate_cell_cycle(cds,
                           g1s_markers = cc.genes$s.genes,
                           g2m_markers = cc.genes$g2m.genes)

# Set Proliferation Index cutoffs
colData(cds)$proliferation_cutoff =
  ifelse(colData(cds)$cell_type == "MCF7", 2.25, 2)

colData(cds)$proliferation_high =
  colData(cds) %>%
  as.data.frame() %>%
  group_by(cell_type) %>%
  mutate(proliferation_high = proliferation_index > proliferation_cutoff) %>%
  pull(proliferation_high)

colData(cds) %>%
  as.data.frame() %>%
  group_by(cell_type, proliferation_high) %>%
  summarise(n = n())


ggplot() +
  geom_density(
    data = colData(cds) %>%
      as.data.frame,
    aes(x = proliferation_index),
    size = .5,
    color = "red"
  ) +
  geom_vline(
    data = colData(cds) %>%
      as.data.frame %>%
      dplyr::select(cell_type, proliferation_cutoff) %>%
      distinct(),
    aes(xintercept = proliferation_cutoff),
    color = "black"
  ) +
  monocle_theme_opts() +
  scale_fill_viridis_c() +
  xlab("Proliferation Index") +
  ylab("Density") +
  facet_wrap( ~ cell_type, ncol = 1) +
  theme(legend.position = "none",
        text = element_text(size = 6)) +
  ggsave(
    "supplemental_figure11_proliferation_ridges.pdf",
    height = 4,
    width = 1.5
  )


test_results =
  colData(cds) %>%
  as.data.frame() %>%
  filter(vehicle) %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(test_res_vehicle = furrr::future_map(data, function(coldata_subset, cds_all) {
    cells_to_test = coldata_subset$cell
    cds_to_test = cds_all[, cells_to_test]
    test_result_subset = fit_models(cds = cds_to_test,
                                    model_formula_str = "~proliferation_high",
                                    verbose = T)
    test_result_coefs = coefficient_table(test_result_subset)
    test_result_coefs
  }, cds, .progress = T)) %>%
  dplyr::select(-data) %>%
  unnest()


# Correct for multiple testing
test_results =
  test_results %>%
  group_by(term) %>%
  mutate(q_value = p.adjust(p_value))

write.table(
  test_results,
  "vehicle_high_vs_low_adjusted_model_test_res.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)
