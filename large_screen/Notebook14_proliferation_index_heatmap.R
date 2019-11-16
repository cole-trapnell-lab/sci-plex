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
  library(viridis)
  library(pheatmap)
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
                "test_drug_dependence.R",
                sep = ""))
  })


# Set DelayedArray Parameters
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size = 1000e7)

# User defined functions

setwd(paste0(path_to_github,
             "large_screen",
             sep = ""))

# Load processed sciPlex CDS's
sciPlex_cds.list = list()
sciPlex_cds.list[["A549"]] = readRDS("A549_24hrs.RDS")
sciPlex_cds.list[["K562"]] = readRDS("K562_24hrs.RDS")
sciPlex_cds.list[["MCF7"]] = readRDS("MCF7_24hrs.RDS")

# Pull in drug annotations to annotate DEG results
drug_annotations = read.table("bin/drugProperties_short.tsv",
                              header = T, sep = "\t")

drug_annotations =
  drug_annotations %>%
  dplyr::rename(catalog_number = Catalog.Number,
                product_name = Product.Name,
                pathway = Pathway,
                target = Target)


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

drug_annotations = left_join(drug_annotations,pathway_annotations, by = "catalog_number" )


for(cell_line in names(sciPlex_cds.list)){

  # Detect the genes that are expressed in the dataset
  sciPlex_cds.list[[cell_line]] <- detect_genes(sciPlex_cds.list[[cell_line]], min_expr = 0.1)

  # Estimate Size Factors for every cell
  sciPlex_cds.list[[cell_line]] <- estimate_size_factors(sciPlex_cds.list[[cell_line]])

  # Estimate Cell Cycle based on marker gene expression for every cell
  sciPlex_cds.list[[cell_line]] <- estimate_cell_cycle(sciPlex_cds.list[[cell_line]],
                                                       g1s_markers = cc.genes$s.genes,
                                                       g2m_markers = cc.genes$g2m.genes)

}

prolifetation.matrix = list()
for(cell_line in names(sciPlex_cds.list)){
  prolifetation.matrix[[cell_line]] =
    colData(sciPlex_cds.list[[cell_line]]) %>%
    as.data.frame() %>%
    group_by(catalog_number,dose_character) %>%
    mutate(median_proliferation_index = median(proliferation_index)) %>%
    ungroup() %>%
    dplyr::select(product_name, dose_character, median_proliferation_index) %>%
    distinct() %>%
    spread(key = dose_character, value = median_proliferation_index ) %>%
    as.data.frame()

  rownames(prolifetation.matrix[[cell_line]]) = prolifetation.matrix[[cell_line]]$product_name

  prolifetation.matrix[[cell_line]] = prolifetation.matrix[[cell_line]][,2:ncol(prolifetation.matrix[[cell_line]])]
  prolifetation.matrix[[cell_line]][,"0"] = prolifetation.matrix[[cell_line]]["Vehicle","0"]

  prolifetation.matrix[[cell_line]] = prolifetation.matrix[[cell_line]][rownames(prolifetation.matrix[[cell_line]]) != "Vehicle",]

  prolifetation.matrix[[cell_line]] = prolifetation.matrix[[cell_line]] %>% as.matrix()

}


ann_colors = list("Pathway" = c("Antioxidant"="aquamarine3",
                                "Apoptotic regulation"="lemonchiffon4",
                                "Cell cycle regulation"="deepskyblue2",
                                "DNA damage & DNA repair"="slategray4",
                                "Epigenetic regulation"="navy",
                                "Focal adhesion signaling"="brown4",
                                "HIF signaling"="darkgreen",
                                "JAK/STAT signaling"="orangered2",
                                "Metabolic regulation" = "gold2",
                                "Neuronal signaling" = "darkolivegreen",
                                "Nuclear receptor signaling" = "chartreuse3",
                                "PKC signaling" = "plum4",
                                "Protein folding & Protein degradation" = "darkorchid2",
                                "TGF/BMP signaling" = "darkcyan",
                                "Tyrosine kinase signaling" = "firebrick1",
                                "Other" ="darkgoldenrod4"))


row_order =
  drug_annotations %>%
  filter(product_name %in% rownames(prolifetation.matrix[["A549"]])) %>%
  arrange(pathway_level_1,pathway_level_2) %>%
  pull(product_name) %>%
  as.character()

all_proliferation_index = rbind(prolifetation.matrix[["A549"]], prolifetation.matrix[["MCF7"]], prolifetation.matrix[["K562"]])
mat_breaks = seq(min(all_proliferation_index),max(all_proliferation_index),length.out =25)

pheatmap(prolifetation.matrix[["A549"]][row_order,],
         file = "supplemental_figure12_A549_proliferation_index_heatmap.pdf",
         annotation_colors = ann_colors,
         annotation_row =
              data.frame(row.names =
                            drug_annotations %>%
                            filter(product_name %in% rownames(prolifetation.matrix[["A549"]])) %>%
                            pull(product_name) %>%
                            as.character(),
                        Pathway =
                            drug_annotations %>%
                            filter(product_name %in% rownames(prolifetation.matrix[["A549"]])) %>%
                             pull(pathway_level_1) %>%
                            as.character()),
         width = 3,
         height = 12,
         breaks = mat_breaks,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 4,
         cellwidth = 2,
         cellheight = 4,
         show_colnames = F,
         show_rownames = F,
         legend = F,
         annotation_legend = F,
         color = viridis(25))



pheatmap(prolifetation.matrix[["MCF7"]][row_order,],
         file = "supplemental_figure12_MCF7_proliferation_index_heatmap.pdf",
         width = 3,
         height = 12,
         breaks = mat_breaks,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 4,
         cellwidth = 2,
         cellheight = 4,
         show_colnames = F,
         show_rownames = F,
         legend = F,
         annotation_legend = F,
         color = viridis(25))


pheatmap(prolifetation.matrix[["K562"]][row_order,],
         file = "supplemental_figure12_K562_proliferation_index_heatmap.pdf",
         width = 3,
         height = 12,
         breaks = mat_breaks,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 4,
         cellwidth = 2,
         cellheight = 4,
         show_colnames = F, 
         show_rownames = F,
         legend = F,
         color = viridis(25))


fraction_in_low.matrix = list()

prolifetation_cutoff = list()

prolifetation_cutoff[["A549"]] = 2
prolifetation_cutoff[["K562"]] = 2
prolifetation_cutoff[["MCF7"]] = 2.25


for(cell_line in names(sciPlex_cds.list)){
  
  colData(sciPlex_cds.list[[cell_line]])$low_proliferation =
    colData(sciPlex_cds.list[[cell_line]])$proliferation_index <= prolifetation_cutoff[[cell_line]]

  fraction_in_low.matrix[[cell_line]] =
  colData(sciPlex_cds.list[[cell_line]]) %>%
    as.data.frame() %>%
    group_by(catalog_number,dose_character) %>%
    add_tally() %>%
    mutate(fraction_in_low = sum(low_proliferation)/n) %>%
    ungroup() %>%
    dplyr::select(fraction_in_low, product_name, dose_character) %>%
    distinct() %>%
    spread(key = dose_character, value = fraction_in_low, fill = 0 ) %>%
    as.data.frame()
  
  rownames(fraction_in_low.matrix[[cell_line]]) = fraction_in_low.matrix[[cell_line]]$product_name
  
  fraction_in_low.matrix[[cell_line]] = fraction_in_low.matrix[[cell_line]][,2:ncol(fraction_in_low.matrix[[cell_line]])]
  fraction_in_low.matrix[[cell_line]][,"0"] = fraction_in_low.matrix[[cell_line]]["Vehicle","0"]
  
  fraction_in_low.matrix[[cell_line]] = fraction_in_low.matrix[[cell_line]][rownames(fraction_in_low.matrix[[cell_line]]) != "Vehicle",]
  
  fraction_in_low.matrix[[cell_line]] = fraction_in_low.matrix[[cell_line]] %>% as.matrix()
  
}


row_order =
  drug_annotations %>%
  filter(product_name %in% rownames(fraction_in_low.matrix[["A549"]])) %>%
  arrange(pathway_level_1,pathway_level_2) %>%
  pull(product_name) %>%
  as.character()


all_fraction_in_low = rbind(fraction_in_low.matrix[["A549"]], fraction_in_low.matrix[["MCF7"]], fraction_in_low.matrix[["K562"]])
mat_breaks = seq(0,1,length.out =25)

colfunc <- colorRampPalette(c("white", "red"))

pheatmap(fraction_in_low.matrix[["A549"]][row_order,],
         file = "supplemental_figure12_A549_frac_low.pdf",
         width = 3,
         height = 12,
         breaks = mat_breaks,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 4,
         cellwidth = 2,
         cellheight = 4,
         show_colnames = F,
         show_rownames = F,
         legend = F,
         annotation_legend = F,         
         color = colfunc(25))



pheatmap(fraction_in_low.matrix[["MCF7"]][row_order,],
         file = "supplemental_figure12_MCF7_frac_low.pdf",
         width = 3,
         height = 12,
         breaks = mat_breaks,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 4,
         cellwidth = 2,
         cellheight = 4,
         show_colnames = F,
         show_rownames = F,
         legend = F,
         annotation_legend = F,
         color = colfunc(25))


pheatmap(fraction_in_low.matrix[["K562"]][row_order,],
         file = "K562_frac_low.pdf",
         width = 3,
         height = 12,
         breaks = mat_breaks,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 4,
         cellwidth = 2,
         cellheight = 4,
         show_colnames = F, 
         show_rownames = F,
         legend = F,
         color = colfunc(25))


low_bin_a549 = 
  colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  group_by(catalog_number,dose_character) %>%
  add_tally() %>%
  mutate(fraction_in_low = sum(low_proliferation)/n) %>%
  ungroup() %>%
  dplyr::select(fraction_in_low, viability, product_name, dose_character, cell_type,proliferation_index) %>%
  distinct() %>%
  mutate(low_bin = cut(fraction_in_low, breaks= seq(0, 1, length.out = 10), labels=F))
    
low_bin_mcf7 = 
  colData(sciPlex_cds.list[["MCF7"]]) %>%
  as.data.frame() %>%
  group_by(catalog_number,dose_character) %>%
  add_tally() %>%
  mutate(fraction_in_low = sum(low_proliferation)/n) %>%
  ungroup() %>%
  dplyr::select(fraction_in_low, viability, product_name, dose_character, cell_type,proliferation_index) %>%
  distinct() %>%
  mutate(low_bin = cut(fraction_in_low, breaks= seq(0, 1, length.out = 10), labels=F))


low_bin_k562 = 
  colData(sciPlex_cds.list[["K562"]]) %>%
  as.data.frame() %>%
  group_by(catalog_number,dose_character) %>%
  add_tally() %>%
  mutate(fraction_in_low = sum(low_proliferation)/n) %>%
  ungroup() %>%
  dplyr::select(fraction_in_low, viability, product_name, dose_character, cell_type,proliferation_index) %>%
  distinct() %>%
  mutate(low_bin = cut(fraction_in_low, breaks= seq(0, 1, length.out = 10), labels=F))

low_bins = 
  rbind(low_bin_a549,low_bin_mcf7,low_bin_k562) %>% 
  filter(!is.na(low_bin))

low_bins$proliferation_cutoff = 
  ifelse(low_bins$cell_type == "MCF7", 2.25, 2)

ggplot() + 
  geom_density(data =low_bins %>%
                 filter(product_name != "Vehicle"),
               aes(x = proliferation_index), fill = "deepskyblue2") +
  geom_density(data =low_bins %>%
                 filter(product_name == "Vehicle"),
               aes(x = proliferation_index), size = .5, color = "red" ) +
    geom_vline(data =low_bins %>%
                 dplyr::select(cell_type, proliferation_cutoff) %>%
                 distinct(),
               aes(xintercept = proliferation_cutoff), color = "black") +
  monocle_theme_opts() + 
  scale_fill_viridis_c() +
  xlab("Proliferation Index") +
  ylab("Density") +
  facet_wrap(~cell_type, ncol = 1)+
  theme(legend.position = "none", 
        text = element_text(size = 6)) +  
  ggsave("proliferation_ridges.pdf", height = 4, width = 1.5)




ggplot(low_bins %>% 
         dplyr::select(-proliferation_index) %>%
         distinct(),
       aes(y = viability, x = as.factor(low_bin))) + 
  geom_boxplot(fill = "grey90", outlier.size = 0, outlier.stroke =  0) +
  geom_jitter(stroke = 0, size =1)+
  geom_point(data = 
               low_bins %>% 
               dplyr::select(-proliferation_index) %>%
               distinct() %>%
               filter(product_name == "Vehicle"),
             color = "red",
             aes(y = viability, x = low_bin)) +
  facet_wrap(~cell_type, ncol = 1) +
  monocle_theme_opts() + 
  xlab("Percent Low Proliferation") +
  ylab("Viability") +
  theme(legend.position = "none", 
        text = element_text(size = 6)) +  
  ggsave("viability_vs_proliferation_boxplot.pdf", height = 5, width = 2.5)



# proliferation_index UMAP
## Chisq_Cutoff_Residual_model_umi_replicate
colData(sciPlex_cds.list[["K562"]]) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = proliferation_index),size = 0.1, stroke = 0) +
  scale_color_viridis_c() +
  xlab("Component 1") +
  ylab("Component 2") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("supplemental_figure12_sciPlex_K562_24hr_UMAP_proliferation_index.png",
         width = 1.25, height = 1.25, unit = "in", dpi = 600)

## Chisq_Cutoff_Residual_model_umi_replicate
colData(sciPlex_cds.list[["A549"]]) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = proliferation_index),size = 0.1, stroke = 0) +
  scale_color_viridis_c() +
  xlab("Component 1") +
  ylab("Component 2") +
  theme_void() +
  theme(legend.position = "none") +
  
  ggsave("supplemental_figure12_sciPlex_A549_24hr_UMAP_proliferation_index.png",
         width = 1.25, height = 1.25, unit = "in", dpi = 600)

## Chisq_Cutoff_Residual_model_umi_replicate
colData(sciPlex_cds.list[["MCF7"]]) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = proliferation_index),size = 0.1, stroke = 0) +
  scale_color_viridis_c() +
  xlab("Component 1") +
  ylab("Component 2") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("supplemental_figure12_sciPlex_MCF7_24hr_UMAP_proliferation_index.png",
         width = 1.25, height = 1.25, unit = "in", dpi = 600)



sciPlex_counts.list <- list()
sciPlex_models.list <- list()
sciPlex_model_stats.list <- list()

for(cell_line in names(sciPlex_cds.list)){
  # Count the number of cells from each well and normalize per replicate and culture plate
  sciPlex_counts.list[[cell_line]] =
    count_cells(sciPlex_cds.list[[cell_line]], normalization_terms=c("cell_type", "culture_plate","replicate"))
  
  sciPlex_models.list[[cell_line]] =
    fit_dose_response_models(sciPlex_counts.list[[cell_line]])
  
  sciPlex_model_stats.list[[cell_line]] =
    extract_dose_model_stats(sciPlex_counts.list[[cell_line]],sciPlex_models.list[[cell_line]] )
  
}

for(cell_line in names(sciPlex_cds.list)){

  # Perform viability estimation on each cell line with the fit sciPlex models
  sciPlex_cds.list[[cell_line]] =
    estimate_viability(sciPlex_cds.list[[cell_line]], sciPlex_models.list[[cell_line]])

  # Set all viability that are estimated as NA to 1
  colData(sciPlex_cds.list[[cell_line]])$viability[is.na(pData(sciPlex_cds.list[[cell_line]])$viability)] = 1
}

viability.compound = list()
for(cell_line in names(sciPlex_cds.list)){
  viability.compound[[cell_line]] =
    colData(sciPlex_cds.list[[cell_line]]) %>%
    as.data.frame() %>%
    dplyr::select(product_name, dose_character, viability) %>%
    distinct() %>%
    pull(product_name)
}

all_viability_compounds =
  union(viability.compound[["A549"]], union(viability.compound[["MCF7"]],viability.compound[["K562"]]))

viability.matrix = list()
for(cell_line in names(sciPlex_cds.list)){
  viability.matrix[[cell_line]] =
    colData(sciPlex_cds.list[[cell_line]]) %>%
    as.data.frame() %>%
    dplyr::select(product_name, dose_character, viability) %>%
    filter(product_name %in% all_viability_compounds) %>%
    distinct() %>%
    spread(key = dose_character, value = viability, fill = 1 ) %>%
    as.data.frame()

  rownames(viability.matrix[[cell_line]]) = viability.matrix[[cell_line]]$product_name

  viability.matrix[[cell_line]] = viability.matrix[[cell_line]][,2:ncol(viability.matrix[[cell_line]])]
  viability.matrix[[cell_line]][,"0"] = viability.matrix[[cell_line]]["Vehicle","0"]

  viability.matrix[[cell_line]] = viability.matrix[[cell_line]][rownames(viability.matrix[[cell_line]]) != "Vehicle",]

  viability.matrix[[cell_line]] = viability.matrix[[cell_line]] %>% as.matrix()

}



ann_colors = list("Pathway" = c("Antioxidant"="aquamarine3",
                                "Apoptotic regulation"="lemonchiffon4",
                                "Cell cycle regulation"="deepskyblue2",
                                "DNA damage & DNA repair"="slategray4",
                                "Epigenetic regulation"="navy",
                                "Focal adhesion signaling"="brown4",
                                "HIF signaling"="darkgreen",
                                "JAK/STAT signaling"="orangered2",
                                "Metabolic regulation" = "gold2",
                                "Neuronal signaling" = "darkolivegreen",
                                "Nuclear receptor signaling" = "chartreuse3",
                                "PKC signaling" = "plum4",
                                "Protein folding & Protein degradation" = "darkorchid2",
                                "TGF/BMP signaling" = "darkcyan",
                                "Tyrosine kinase signaling" = "firebrick1",
                                "Other" ="darkgoldenrod4"))


row_order =
  drug_annotations %>%
  filter(product_name %in% rownames(viability.matrix[["A549"]])) %>%
  arrange(pathway_level_1,pathway_level_2) %>%
  pull(product_name) %>%
  as.character()

all_viability_matrix = rbind(viability.matrix[["A549"]], viability.matrix[["MCF7"]], viability.matrix[["K562"]])

num_breaks = 25
mat_breaks = seq(min(all_viability_matrix),max(all_viability_matrix),length.out =num_breaks)


pheatmap(viability.matrix[["A549"]][row_order,],
         file = "supplemental_figure12_A549_viability_index.pdf",
         annotation_colors = ann_colors,
         width = 3,
         height = 12,
         breaks = mat_breaks,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 4,
         cellwidth = 2,
         cellheight = 4,
         show_colnames = F,
         show_rownames = F,
         legend = F,
         annotation_legend = F,
         color = viridis(option = "magma",25))



pheatmap(viability.matrix[["MCF7"]][row_order,],
         file = "supplemental_figure12_MCF7_viability_index.pdf",
         width = 3,
         height = 12,
         breaks = mat_breaks,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 4,
         cellwidth = 2,
         cellheight = 4,
         show_colnames = F,
         show_rownames = F,
         legend = F,
         annotation_legend = F,
         color = viridis(option = "magma",25))


pheatmap(viability.matrix[["K562"]][row_order,],
         file = "supplemental_figure12_K562_viability_index.pdf",
         width = 3,
         height = 12,
         breaks = mat_breaks,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 4,
         cellwidth = 2,
         cellheight = 4,
         show_colnames = F,         
         legend=F,
         color = viridis(option = "magma",25))
