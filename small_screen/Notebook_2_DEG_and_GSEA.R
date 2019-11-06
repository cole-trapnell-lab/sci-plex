# Input path to dropbox
path_to_github = "~/sci-Plex/"
path_to_monocle3 = paste0(path_to_github,
                          "monocle3_d4a9a35/monocle3/",
                          sep = "")

suppressPackageStartupMessages({
  library(tidymodels)
  library(devtools)
  load_all(path_to_monocle3)
  library(furrr)
  library(viridis)
  library(piano)
  library(UpSetR)
  library(snowfall)
  plan(multicore)
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
  source(paste0(bin_directory,
                "GSA_helper_functions.R",
                sep = ""))
  source(paste0(bin_directory,
                "loadGSCSafe.R",
                sep = ""))
  
})


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
sciPlex_cds = estimate_cell_cycle(sciPlex_cds, cc.genes$s.gene, cc.genes$g2m.gene)


colData(sciPlex_cds)$log_hash_umis = log(colData(sciPlex_cds)$hash_umis)


sciPlex_cds = detect_genes(sciPlex_cds, min_expr = 0.1)
sciPlex_cds = estimate_size_factors(sciPlex_cds)

expressed_genes = row.names(fData(sciPlex_cds)[Matrix::rowSums(as.matrix(exprs(sciPlex_cds)) > 0) > 50 , ])

########## DEG analysis:


vehicle_cells = pData(sciPlex_cds) %>% filter(vehicle == TRUE)
vehicle_cells = as.character(vehicle_cells$cell)

reference_cells = sample_n(pData(sciPlex_cds), 1000)$cell

test_drug_dependence <- function(cds_subset, cds, reference_cells=NULL, min_fraction_cells=0.05, pseudodose=0.01, residualModelFormulaTerms=NULL, cores=1){
    #print(cell_ids_to_test)
    cell_ids_to_test = as.character(cds_subset$cell)

    print(length(cell_ids_to_test))

    cds_subset = cds[,reference_cells]
    cds_subset = detect_genes(cds_subset)
    genes_in_reference_set = row.names(subset(fData(cds_subset), num_cells_expressed > ncol(cds_subset)*min_fraction_cells))
    
    cds_subset = cds[,cell_ids_to_test]
    cds_subset = detect_genes(cds_subset)
    genes_in_treated_set = row.names(subset(fData(cds_subset), num_cells_expressed > ncol(cds_subset)*min_fraction_cells))
    print(length(genes_in_treated_set))
    cds_subset = cds[,base::union(cell_ids_to_test, reference_cells)]
    cds_subset@expression_family = "quasipoisson"

    cds_subset = cds_subset[base::union(genes_in_reference_set, genes_in_treated_set),]
    tmnt = unique(pData(cds_subset)$treatment)[1]
    ct = unique(pData(cds_subset)$cell_type)
    cell_treated = as.character(pData(cds_subset)$treatment) == as.character(tmnt)
    dose_of_drug = log((pData(cds_subset)$dose * cell_treated)  + pseudodose)
    drug_var_name = paste("dose_", tmnt, sep="")
    pData(cds_subset)[,drug_var_name] = dose_of_drug
    message(paste(ncol(cds_subset), ct, "cells treated with", tmnt))
    message(paste("Fitting drug/dose models to ", nrow(cds_subset), "genes"))
    modelFormulaStr = paste("~",  paste("dose_",tmnt, sep=""))
    
    if (is.null(residualModelFormulaTerms) == FALSE & length(residualModelFormulaTerms) > 0){
      for (i in (1:length(residualModelFormulaTerms))){
        if (length(unique(pData(cds_subset)[,residualModelFormulaTerms[i]])) > 1){
          modelFormulaStr = paste(modelFormulaStr, "+", residualModelFormulaTerms[i])
        }
      }
    }
    if (cds_subset@expression_family %in% c("zipoisson", "zinegbinomial"))
      modelFormulaStr = paste(modelFormulaStr, "| log_hash_umis")
    print (modelFormulaStr)

    sci_plex_model_tbl = fit_models(cds_subset, model_formula_str = modelFormulaStr, cores=cores, reltol=1e-5, verbose=TRUE)
    sci_plex_model_coefs = coefficient_table(sci_plex_model_tbl)
    sci_plex_model_coefs
}

plan("multiprocess")
old_omp_num_threads = Sys.getenv("OMP_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
simple_model_test_res =
  colData(small_cds) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  mutate(
    test_res = purrr::map(
      data,
      .f = test_drug_dependence,
      small_cds,
      union(vehicle_cells, reference_cells),
      residualModelFormulaTerms = NULL,
      cores = 1
    )
  )

simple_model_test_res = 
  simple_model_test_res %>% 
  dplyr::select(-data) %>% 
  unnest()

Sys.setenv(OMP_NUM_THREADS = old_omp_num_threads)


# Correct for multiple testing
simple_model_test_res =
  simple_model_test_res %>%
  group_by(term) %>%
  mutate(q_value = p.adjust(p_value))

# Can be downloaded from website found under Supplemental Table 2
## Insert link here
write.table(simple_model_test_res, 
            "A549_small_screen_DEGs.txt", 
            row.names=F, 
            quote=F, 
            sep="\t")



simple_model_hit_table = simple_model_test_res %>%
  filter(grepl("dose_", term)) %>%
  mutate(hit = ifelse(q_value < 0.05, 1, 0)) %>%
  dplyr::select(id, term, hit) %>%
  spread(term, hit, fill = 0)


colnames(simple_model_hit_table) = stringr::str_replace(colnames(adjusted_model_hit_table), "dose_", "")
pdf(file = "supplemental_figure3_hits_upset.pdf",
    onefile = FALSE,
    height = 3) # or other device
upset(
  as.data.frame(simple_model_hit_table),
  mainbar.y.label = "DEGs in intersection",
  sets.x.label = "Total DEGs"
)
dev.off()

########## GSEA ##########

## Load Gene Set Collection
hallmarks <-loadGSCSafe(file= paste0(path_to_github,
                                     "GeneSets/h.all.v6.0.symbols.gmt",
                                     sep = ""))

# This function is designed to collect GSA stats on each drug + cell_type combo
drug_effect_gsa <- function(drug_gene_test_res, gsc, min_genes_in_set=5, ncores=1){
    drug_gene_test_res = drug_gene_test_res %>% 
             group_by(gene_short_name) %>%
             filter(abs(test_val) == max(abs(test_val))) 
    directions=data.frame(row.names=as.character(drug_gene_test_res$gene_short_name), stats=drug_gene_test_res$test_val)
    geneLevelStats=data.frame(row.names=as.character(drug_gene_test_res$gene_short_name), stats=drug_gene_test_res$test_val)
    tryCatch({
        gsares <- runGSA(geneLevelStats=geneLevelStats, direction=geneLevelStats, geneSetStat="mean", gsc=gsc, nPerm=10000, ncpu=ncores, gsSizeLim=c(min_genes_in_set, Inf), verbose=T)
        return (gsares)
    }, error = function(e) { print(e); return (NA) } )
}
Sys.setenv(OMP_NUM_THREADS = 1)

simple_model_gsa_res = simple_model_test_res %>% filter(grepl("dose_", term)) %>% group_by(treatment, term) %>% 
  nest(-term) %>% 
  mutate(
    gsares = purrr::map(.f = drug_effect_gsa, .x = data, hallmarks, ncores=4)
  ) 

# Extract p values and GSA test statistics into a table
simple_model_gsa_summary = 
  simple_model_gsa_res %>%
      mutate(gsa_summary = purrr::map(.f = GSAsummaryTable, .x = gsares)) %>% 
      unnest(gsa_summary)


## Make a heatmap showing the effects of all the drugs on one cell line:
simple_model_gsa_summary = simple_model_gsa_summary %>%
  mutate(up_stat = `Stat (mix.dir.up)`, dn_stat = `Stat (mix.dir.dn)`) %>%
  dplyr::select(term,
                Name,
                `Genes (tot)`,
                `p (mix.dir.up)`,
                up_stat,
                `p (mix.dir.dn)`,
                dn_stat)
simple_model_gsa_summary = simple_model_gsa_summary %>% mutate(gsa_stat = ifelse((is.na(up_stat) == FALSE &
                                                                                    `p (mix.dir.up)` < `p (mix.dir.dn)`),
                                                                                 up_stat,
                                                                                 ifelse(is.na(dn_stat) == FALSE,-dn_stat, 0)
))

simple_model_gsa_stat = simple_model_gsa_summary %>% dplyr::select(Name, term, gsa_stat) %>% spread(term, gsa_stat)

gsa_effects = simple_model_gsa_stat

# Get the top gene sets for both
num_top_gene_sets = 4

top_simple_gene_sets =
  simple_model_gsa_summary %>%
  group_by(term) %>%
  top_n(num_top_gene_sets, abs(gsa_stat))

top_simple_gene_sets =
  top_simple_gene_sets %>%
  arrange(term, gsa_stat) %>%
  dplyr::select(Name, term, gsa_stat) %>%
  mutate(Name = stringr::str_split_fixed(Name, "%", 3)[, 1]) %>% as.data.frame
sig_pathway_names = unique(as.character(top_simple_gene_sets$Name))


# Make into a matrix
gsa_effect_matrix = gsa_effects %>% dplyr::select(-Name) %>% as.matrix
gsa_effect_matrix[is.na(gsa_effect_matrix)] = 0
row.names(gsa_effect_matrix) = stringr::str_split_fixed(gsa_effects$Name, "%", 3)[, 1]

# Make into a data frame
gsa_effect_df =
  gsa_effects %>%
  gather(key = term, value = gsa_stat, contains("dose"))

gsa_effect_df[is.na(gsa_effect_df)] = 0
gsa_effect_df$Name = stringr::str_split_fixed(gsa_effects$Name, "HALLMARK_", 2)[, 2]
gsa_effect_df$Name = stringr::str_replace(gsa_effect_df$Name, "_", " ")


gsa_effect_matrix[gsa_effect_matrix < -10] = -10
gsa_effect_matrix[gsa_effect_matrix > 10] = 10

matrix_to_plot = gsa_effect_matrix[sig_pathway_names, ]

colnames(matrix_to_plot)  = (colnames(matrix_to_plot) %>% stringr::str_split_fixed("_", 2))[, 2]

rn =   stringr::str_split_fixed(rownames(matrix_to_plot), pattern = "HALLMARK_", n =  2)[, 2]
rownames(matrix_to_plot) = stringr::str_replace_all(rn, "_", " ")
pheatmap(
  matrix_to_plot,
  cex = 0.5,
  filename = "supplemental_figure3_sciPlex_GOheatmap.pdf",
  width = 3,
  height = 2.25,
  cellwidth = 8,
  cellheight = 8,
  fontsize = 10,
  fontsize_row = 25,
  fontsize_col = 25,
  legend = F,
  treeheight_row = 15,
  treeheight_col = 15,
  border_color = "black"
)


