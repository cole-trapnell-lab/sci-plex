# devtools::install_github(repo = 'cole-trapnell-lab/monocle3', ref = "d4a9a35")

# Input path to github directory
path_to_github = "~/sci-Plex/"

# Set directory for sciPlex bin
bin_directory = paste0(path_to_github,
                       "bin",
                       sep = "")

suppressPackageStartupMessages({
  library(ggplot2)
  library(scran)
  library(piano)
  library(dplyr)
  library(tidymodels)
  library(snowfall)
  library(viridis)
  library(pheatmap)
  library(monocle3)
})

## User defined function
append_umap_coordinates = function(cds){
  num_dim = dim(cds@reducedDims[["UMAP"]])[2]
  for (i in seq(1,num_dim)){
    new_col_name = paste("umap",i,sep = "")
    colData(cds)[[new_col_name]] = cds@reducedDims[["UMAP"]][,i]
  }
  return(cds)
}

## Return an effects matrix from DEG test result table
compute_effect_matrix = function(test_res){
  dose_terms =
    test_res %>%
    dplyr::filter(grepl("dose_", term))

  effect_matrix =
    dose_terms %>%
    ungroup() %>%
    dplyr::select(id, treatment, estimate) %>%
    spread(treatment, estimate, fill=0)

  return(effect_matrix)
}


# Load previously computed DEG results
A549_test_res = read.table("Supplementary_Table_5_A549.txt",
                              header = T)

K562_test_res = read.table("Supplementary_Table_5_K562.txt",
                              header = T)

MCF7_test_res = read.table("Supplementary_Table_5_MCF7.txt",
                              header = T)


drug_info_pdata = data.frame()

A549_effect_matrix = compute_effect_matrix(A549_test_res)
colnames(A549_effect_matrix) = c("id", paste("A549", colnames(A549_effect_matrix)[-1], sep="_"))

MCF7_effect_matrix = compute_effect_matrix(MCF7_test_res)
colnames(MCF7_effect_matrix) = c("id", paste("MCF7", colnames(MCF7_effect_matrix)[-1], sep="_"))

K562_effect_matrix = compute_effect_matrix(K562_test_res)
colnames(K562_effect_matrix) = c("id", paste("K562", colnames(K562_effect_matrix)[-1], sep="_"))

effect_matrix = full_join(A549_effect_matrix, full_join(MCF7_effect_matrix, K562_effect_matrix, by="id", fill=0), fill=0, by="id")

gene_ids = as.character(effect_matrix$id)
effect_matrix = effect_matrix[,-1]
effect_matrix = as.matrix(effect_matrix)
effect_matrix[is.na(effect_matrix)] = 0
row.names(effect_matrix) = gene_ids

# Pull in drug annotations to annotate DEG results
drug_annotations = read.table("bin/drugProperties_short.tsv",
                              header = T, sep = "\t")

drug_annotations =
  drug_annotations %>%
  dplyr::rename(catalog_number = Catalog.Number,
                product_name = Product.Name,
                pathway = Pathway,
                target = Target)


# Pull in drug pathway annotations to annotate results
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

row.names(drug_annotations) = drug_annotations$catalog_number



cell_type_drug_names = stringr::str_split_fixed(colnames(effect_matrix), "_", 2)
drug_info_pdata = drug_annotations[cell_type_drug_names[,2],]
drug_info_pdata$cell_type = cell_type_drug_names[,1]

A549_dose_terms = A549_test_res %>% filter(grepl("dose_", term))
MCF7_dose_terms = MCF7_test_res %>% filter(grepl("dose_", term))
K562_dose_terms = K562_test_res %>% filter(grepl("dose_", term))

deg_summary =
  bind_rows(A549_dose_terms, MCF7_dose_terms, K562_dose_terms) %>%
  filter(q_value < 0.1) %>%
  group_by(cell_type, treatment) %>%
  summarize(num_degs = n())

drug_info_pdata = left_join(drug_info_pdata, deg_summary, by=c("cell_type"="cell_type", "catalog_number"="treatment"))

colnames(effect_matrix) = paste0(drug_info_pdata$cell_type,"_",drug_info_pdata$catalog_number)
rownames(drug_info_pdata) = paste0(drug_info_pdata$cell_type,"_",drug_info_pdata$catalog_number)


effect_cds = new_cell_data_set(effect_matrix, drug_info_pdata )

colData(effect_cds)$Size_Factor = 1



A549_deg_ids = A549_test_res %>% filter(grepl("dose_", term) & q_value < 0.1)
A549_deg_ids = unique(A549_deg_ids$id)
K562_deg_ids = K562_test_res %>% filter(grepl("dose_", term) & q_value < 0.1)
K562_deg_ids = unique(K562_deg_ids$id)
MCF7_deg_ids = MCF7_test_res %>% filter(grepl("dose_", term) & q_value < 0.1)
MCF7_deg_ids = unique(MCF7_deg_ids$id)

deg_ids = union(union(A549_deg_ids, K562_deg_ids), MCF7_deg_ids)

effect_cds = effect_cds[deg_ids,colData(effect_cds)$product_name != "Vehicle"]

effect_cds = preprocess_cds(effect_cds,
                      method = 'PCA',
                      num_dim = 25,
                      norm_method = 'size_only',
                      verbose = T)

cells_list <- list()
cell_types <- unique(colData(effect_cds)[,"cell_type"])

for(genotype in cell_types){
  cells_list[[genotype]] <- row.names(colData(effect_cds)[colData(effect_cds)[,"cell_type"] == genotype,])
}

PCA_cells_1 = effect_cds@reducedDims[["PCA"]][cells_list[[1]],]
PCA_cells_2 = effect_cds@reducedDims[["PCA"]][cells_list[[2]],]
PCA_cells_3 = effect_cds@reducedDims[["PCA"]][cells_list[[3]],]



corrected_PCA = scran::mnnCorrect(t(PCA_cells_1),
                                    t(PCA_cells_2),
                                    t(PCA_cells_3))

corrected_PCA = t(do.call("cbind", corrected_PCA[["corrected"]]))

row.names(corrected_PCA) = c(row.names(PCA_cells_1),
                              row.names(PCA_cells_2),
                              row.names(PCA_cells_3))

corrected_PCA = corrected_PCA[colnames(effect_cds),]
effect_cds@reducedDims[["PCA"]]= as.matrix(corrected_PCA)



effect_cds = reduce_dimension(effect_cds,
                       max_components = 2,
                       reduction_method = 'UMAP',
                       umap.metric = 'cosine',
                       umap.n_neighbors = 10,
                       umap.min_dist = 0.1,
                       verbose = TRUE)

effect_cds = append_umap_coordinates(effect_cds)

effect_cds = cluster_cells(effect_cds, method="louvain", k=10)

pathway_colors =  c("Antioxidant"="aquamarine3",
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
                    "Other" ="darkgoldenrod4")

colData(effect_cds) %>%
  as.data.frame() %>%
  ggplot(aes(x=umap1, y=umap2, color=pathway_level_1, size=I(0.35), shape=cell_type)) +
  geom_point() +
  theme(legend.position = "right", text = element_text(size = 6),  legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.5,"line")) +
  guides(guides(colour = guide_legend(override.aes = list(size=1)))) +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  scale_color_manual("Pathway", values = pathway_colors) +
  theme(legend.position = "none") +
  ggsave("supplemental_figure19_drug_LM_betas_UMAP.pdf", height = 2, width = 2, units = "in", dpi=600)

#############


coldata_umap_matrix =
  colData(effect_cds) %>%
  as.data.frame() %>%
  dplyr::select(umap1, umap2) %>%
  as.matrix()

coldata_umap_matrix =
  effect_cds@reducedDims[["PCA"]]
  


contrasts = colData(effect_cds) %>%
  as.data.frame() %>%
  dplyr::select(catalog_number, cell_type, umap1, umap2) %>%
  mutate(row_names = rownames(colData(effect_cds))) %>%
  group_by(catalog_number) %>%
  nest() %>%
  mutate(combos =  purrr::map(data, .f = function(pdata_subset){
    df =
    pdata_subset %>%
      tidyr::expand(cell_type,cell_type) %>%
      filter(cell_type > cell_type1)

    tibble(cell_type1 =df$cell_type,
            cell_type2 =df$cell_type1)

  } ))

contrasts =
  contrasts %>%
  dplyr::select(-data) %>%
  unnest()

  contrasts =
  contrasts %>%
  mutate(cell_type1 = paste(cell_type1,catalog_number,  sep = "_"),
        cell_type2 = paste(cell_type2, catalog_number, sep = "_"))

contrasts$umap_distance =
  mapply(function(idx1, idx2){
      dist(rbind(coldata_umap_matrix[idx1,], coldata_umap_matrix[idx2,]))
  }, contrasts$cell_type1, contrasts$cell_type2)

contrasts_matrix =
  contrasts %>%
  mutate( cell_type1 = stringr::str_split_fixed(contrasts$cell_type1,"_",n = 2)[,1],
          cell_type2 = stringr::str_split_fixed(contrasts$cell_type2,"_",n = 2)[,1]) %>%
  mutate(combo = paste(cell_type1, cell_type2)) %>%
  dplyr::select(umap_distance, catalog_number, combo) %>%
  spread(key = combo, value = umap_distance)

catalog_number_contrasts_matrix =
  data.frame(catalog_number = contrasts_matrix$catalog_number)

contrasts_matrix =
  contrasts_matrix[,2:ncol(contrasts_matrix)] %>%
  as.matrix()

catalog_number_contrasts_matrix =
  left_join(catalog_number_contrasts_matrix, drug_annotations, by = "catalog_number")

rownames(contrasts_matrix) = catalog_number_contrasts_matrix$product_name




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
                              "Other" ="darkgoldenrod4",
                              "Vehicle"="grey80"))
row_order =
  drug_annotations %>%
  filter(product_name %in% rownames(contrasts_matrix) )%>%
  arrange(pathway_level_1,pathway_level_2) %>%
  pull(product_name) %>%
  as.character()

contrasts_matrix = contrasts_matrix[row_order,]
mat_breaks = seq(min(contrasts_matrix),sort(contrasts_matrix,decreasing = T)[10],length.out =25)

cuts = 10

ph = pheatmap(contrasts_matrix,
         file = "supplemental_figure20_pairwise_distance_UMAP_on_betas_pca.pdf",
         width = 6,
         height = 12,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = F,
         show_colnames = F,
         breaks = mat_breaks,
         clustering_method = "ward.D2",
         legend = F,
         cellwidth = 5,
         cellheight = 2.55,
         annotation_legend = F,
         fontsize = 6,
         treeheight_row = 0,
         color = viridis(25))


contrasts_matrix_cluster = cutree(ph$tree_row, k = cuts,)

ph = pheatmap(contrasts_matrix,
              file = "supplemental_figure20_pairwise_distance_UMAP_on_betas_cut_rows_PCA.pdf",
              width = 6,
              height = 12,
              cluster_cols = F,
              cluster_rows = T,
              show_rownames = F,
              show_colnames = T,
              breaks = mat_breaks,
              clustering_method = "ward.D2",
              legend = F,
              cellwidth = 4,
              cellheight = 2.5,
              cutree_rows = 10,
              annotation_row =
                data.frame(row.names = names(contrasts_matrix_cluster),
                           cluster = as.factor(contrasts_matrix_cluster)),
              annotation_legend = T,
              fontsize = 8,
              treeheight_row = 0,
              color = viridis(25))


pheatmap(contrasts_matrix[contrasts_matrix_cluster == 1, ],
         file = "supplemental_figure20_pairwise_distance_UMAP_on_betas_same_pca.pdf",
         annotation_colors = ann_colors,
         annotation_row =
              data.frame(row.names =
                            catalog_number_contrasts_matrix %>%
                            filter(product_name %in% rownames(contrasts_matrix)) %>%
                            pull(product_name) %>%
                            as.character(),
                        Pathway =
                            catalog_number_contrasts_matrix %>%
                            filter(product_name %in% rownames(contrasts_matrix)) %>%
                             pull(pathway_level_1) %>%
                            as.character()),
         width = 4,
         height = 12,
         cluster_cols = F,
         cluster_rows = F,
         legend = F,
         show_colnames = F,
         annotation_legend = F,
         annotation_names_row = F,
         breaks = mat_breaks,
         clustering_method = "ward.D2",
         cellwidth = 6,
         cellheight = 6,
         fontsize = 7,
         treeheight_row = 0,
         color = viridis(25))

pheatmap(contrasts_matrix[contrasts_matrix_cluster == 9,],
         file = "supplemental_figure20_pairwise_distance_UMAP_on_betas_K562_pca.pdf",
         annotation_colors = ann_colors,
         annotation_row =
              data.frame(row.names =
                            catalog_number_contrasts_matrix %>%
                            filter(product_name %in% rownames(contrasts_matrix)) %>%
                            pull(product_name) %>%
                            as.character(),
                        Pathway =
                            catalog_number_contrasts_matrix %>%
                            filter(product_name %in% rownames(contrasts_matrix)) %>%
                             pull(pathway_level_1) %>%
                            as.character()),
         width = 4,
         height = 6,
         cluster_cols = F,
         cluster_rows = F,
         legend = F,
         annotation_legend = F,
         annotation_names_row = F,
         show_colnames = F,
         breaks = mat_breaks,
         cellwidth = 6,
         cellheight = 6,
         fontsize = 7,
         treeheight_row = 0,
         color = viridis(25))

pheatmap(contrasts_matrix[contrasts_matrix_cluster == 6,],
         file = "pairwise_distance_UMAP_on_betas_A549_pca.pdf",
         annotation_colors = ann_colors,
         annotation_row =
              data.frame(row.names =
                            catalog_number_contrasts_matrix %>%
                            filter(product_name %in% rownames(contrasts_matrix)) %>%
                            pull(product_name) %>%
                            as.character(),
                        Pathway =
                            catalog_number_contrasts_matrix %>%
                            filter(product_name %in% rownames(contrasts_matrix)) %>%
                             pull(pathway_level_1) %>%
                            as.character()),
         width = 4,
         height = 6,
         cluster_cols = F,
         cluster_rows = F,
         legend = F,
         annotation_names_row = F,
         annotation_legend = F,
         show_colnames = F,
         breaks = mat_breaks,
         cellwidth = 6,
         cellheight = 6,
         fontsize = 7,
         treeheight_row = 0,
         color = viridis(25))


pheatmap(contrasts_matrix[contrasts_matrix_cluster == 9,],
         file = "supplemental_figure20_pairwise_distance_UMAP_on_betas_A549_2.pdf",
         annotation_colors = ann_colors,
         annotation_row =
           data.frame(row.names =
                        catalog_number_contrasts_matrix %>%
                        filter(product_name %in% rownames(contrasts_matrix)) %>%
                        pull(product_name) %>%
                        as.character(),
                      Pathway =
                        catalog_number_contrasts_matrix %>%
                        filter(product_name %in% rownames(contrasts_matrix)) %>%
                        pull(pathway_level_1) %>%
                        as.character()),
         width = 4,
         height = 6,
         cluster_cols = F,
         cluster_rows = F,
         legend = F,
         annotation_names_row = F,
         annotation_legend = F,
         show_colnames = F,
         breaks = mat_breaks,
         cellwidth = 6,
         cellheight = 6,
         fontsize = 7,
         treeheight_row = 0,
         color = viridis(25))

pheatmap(contrasts_matrix[contrasts_matrix_cluster == 5,],
         file = "supplemental_figure20_pairwise_distance_UMAP_on_betas_MCF7_pca.pdf",
         annotation_colors = ann_colors,
         annotation_row =
              data.frame(row.names =
                            catalog_number_contrasts_matrix %>%
                            filter(product_name %in% rownames(contrasts_matrix)) %>%
                            pull(product_name) %>%
                            as.character(),
                        Pathway =
                            catalog_number_contrasts_matrix %>%
                            filter(product_name %in% rownames(contrasts_matrix)) %>%
                             pull(pathway_level_1) %>%
                            as.character()),
         width = 4,
         height = 6,
         cluster_cols = F,
         cluster_rows = F,
         legend = F,
         annotation_legend = F,
         annotation_names_row = F,
         show_colnames = F,
         breaks = mat_breaks,
         cellwidth = 6,
         cellheight = 6,
         fontsize = 7,
         treeheight_row = 0,
         color = viridis(25))



pheatmap(contrasts_matrix[contrasts_matrix_cluster == 10,],
         file = "supplemental_figure_20_pairwise_distance_UMAP_on_betas_MCF7_2.pdf",
         annotation_colors = ann_colors,
         annotation_row =
           data.frame(row.names =
                        catalog_number_contrasts_matrix %>%
                        filter(product_name %in% rownames(contrasts_matrix)) %>%
                        pull(product_name) %>%
                        as.character(),
                      Pathway =
                        catalog_number_contrasts_matrix %>%
                        filter(product_name %in% rownames(contrasts_matrix)) %>%
                        pull(pathway_level_1) %>%
                        as.character()),
         width = 4,
         height = 6,
         cluster_cols = F,
         cluster_rows = F,
         legend = F,
         annotation_legend = F,
         annotation_names_row = F,
         show_colnames = F,
         breaks = mat_breaks,
         cellwidth = 6,
         cellheight = 6,
         fontsize = 7,
         treeheight_row = 0,
         color = viridis(25))


ggplot(data = colData(effect_cds) %>%
         as.data.frame()) +
  geom_point(aes(x=umap1, y=umap2, shape=cell_type), color= "grey80", size=0.35) +
  geom_point(data = colData(effect_cds) %>%
               as.data.frame() %>%
               dplyr::filter(catalog_number == "S1090"),
             aes(x=umap1, y=umap2, shape=cell_type, color=pathway_level_1), size=2,) +
  theme(legend.position = "none") +
  theme_void()+
  scale_color_manual("Pathway", values = pathway_colors) +
  theme(legend.position = "none") +
  ggsave("abexinostat_UMAP.pdf", height = 2, width = 2, units = "in", dpi=600)


ggplot(data = colData(effect_cds) %>%
         as.data.frame()) +
  geom_point(aes(x=umap1, y=umap2, shape=cell_type), color= "grey80", size=0.35) +
  geom_point(data = colData(effect_cds) %>%
               as.data.frame() %>%
               dplyr::filter(catalog_number == "S1085"),
             aes(x=umap1, y=umap2, shape=cell_type, color=pathway_level_1), size=2,) +
  theme(legend.position = "none") +
  theme_void()+
  scale_color_manual("Pathway", values = pathway_colors) +
  theme(legend.position = "none") +
  ggsave("supplemental_figure_20_belinostat_UMAP.pdf", height = 2, width = 2, units = "in", dpi=600)

  ggplot(data = colData(effect_cds) %>%
    as.data.frame()) +
  geom_point(aes(x=umap1, y=umap2, shape=cell_type), color= "grey80", size=0.35) +
  geom_point(data = colData(effect_cds) %>%
                as.data.frame() %>%
                dplyr::filter(catalog_number == "S1033"),
          aes(x=umap1, y=umap2, shape=cell_type, color=pathway_level_1), size=2,) +
  theme(legend.position = "none") +
  theme_void()+
  scale_color_manual("Pathway", values = pathway_colors) +
  theme(legend.position = "none") +
  ggsave("supplemental_figure_20_nilotinib_UMAP.pdf", height = 2, width = 2, units = "in", dpi=600)


ggplot(data = colData(effect_cds) %>%
  as.data.frame()) +
geom_point(aes(x=umap1, y=umap2, shape=cell_type), color= "grey80", size=0.35) +
geom_point(data = colData(effect_cds) %>%
              as.data.frame() %>%
              dplyr::filter(catalog_number == "S1014"),
        aes(x=umap1, y=umap2, shape=cell_type, color=pathway_level_1), size=2,) +
theme(legend.position = "none") +
theme_void()+
scale_color_manual("Pathway", values = pathway_colors) +
theme(legend.position = "none") +
ggsave("supplemental_figure_20_bosutnib_UMAP.pdf", height = 2, width = 2, units = "in", dpi=600)



ggplot(data = colData(effect_cds) %>%
         as.data.frame()) +
  geom_point(aes(x=umap1, y=umap2, shape=cell_type), color= "grey80", size=0.35) +
  geom_point(data = colData(effect_cds) %>%
               as.data.frame() %>%
               dplyr::filter(catalog_number == "S1628"),
             aes(x=umap1, y=umap2, shape=cell_type, color=pathway_level_1), size=2,) +
  theme(legend.position = "none") +
  theme_void()+
  scale_color_manual("Pathway", values = pathway_colors) +
  theme(legend.position = "none") +
  ggsave("supplemental_figure_20_triamcinolone_UMAP.pdf", height = 2, width = 2, units = "in", dpi=600)


ggplot(data = colData(effect_cds) %>%
         as.data.frame()) +
  geom_point(aes(x=umap1, y=umap2, shape=cell_type), color= "grey80", size=0.35) +
  geom_point(data = colData(effect_cds) %>%
               as.data.frame() %>%
               dplyr::filter(catalog_number == "S2736"),
             aes(x=umap1, y=umap2, shape=cell_type, color=pathway_level_1), size=2,) +
  theme(legend.position = "none") +
  theme_void()+
  scale_color_manual("Pathway", values = pathway_colors) +
  theme(legend.position = "none") +
  ggsave("supplemental_figure_20_TG101_UMAP.pdf", height = 2, width = 2, units = "in", dpi=600)



ggplot(data = colData(effect_cds) %>%
         as.data.frame()) +
  geom_point(aes(x=umap1, y=umap2, shape=cell_type), color= "grey80", size=0.35) +
  geom_point(data = colData(effect_cds) %>%
               as.data.frame() %>%
               dplyr::filter(catalog_number == "S2673"),
             aes(x=umap1, y=umap2, shape=cell_type, color=pathway_level_1), size=2,) +
  theme(legend.position = "none") +
  theme_void()+
  scale_color_manual("Pathway", values = pathway_colors) +
  theme(legend.position = "none") +
  ggsave("supplemental_figure_20_trametinib_UMAP.pdf", height = 2, width = 2, units = "in", dpi=600)


ggplot(data = colData(effect_cds) %>%
         as.data.frame()) +
  geom_point(aes(x=umap1, y=umap2, shape=cell_type), color= "grey80", size=0.35) +
  geom_point(data = colData(effect_cds) %>%
               as.data.frame() %>%
               dplyr::filter(catalog_number == "S1191"),
             aes(x=umap1, y=umap2, shape=cell_type, color=pathway_level_1), size=2,) +
  theme(legend.position = "none") +
  theme_void()+
  scale_color_manual("Pathway", values = pathway_colors) +
  theme(legend.position = "none") +
  ggsave("supplemental_figure_20_Fulvestrant_UMAP.pdf", height = 2, width = 2, units = "in", dpi=600)

gene_meta_data = bind_rows(A549_test_res %>% filter(grepl("dose_", term) & q_value < 0.1) %>% ungroup() %>% dplyr::select(id, gene_short_name),
                           K562_test_res %>% filter(grepl("dose_", term) & q_value < 0.1)%>% ungroup() %>% dplyr::select(id, gene_short_name),
                           MCF7_test_res %>% filter(grepl("dose_", term) & q_value < 0.1)%>% ungroup() %>% dplyr::select(id, gene_short_name))
gene_meta_data = gene_meta_data %>% distinct()

plot_cell_clusters(effect_cds, color_by = "louvain_component", x = 1, y = 2, cell_size=0.5)


##############


A549_degs_for_plot =
  A549_test_res %>%
  filter(grepl("dose_", term) & q_value < 0.05) %>%
  dplyr::select(cell_type, term, gene_short_name) %>%
  dplyr::mutate(catalog_number = stringr::str_split_fixed(term, "_", 2)[,2] )

K562_degs_for_plot =
  K562_test_res %>%
  filter(grepl("dose_", term) & q_value < 0.05) %>%
  dplyr::select(cell_type, term, gene_short_name) %>%
  dplyr::mutate(catalog_number = stringr::str_split_fixed(term, "_", 2)[,2] )

MCF7_degs_for_plot =
  MCF7_test_res %>%
  filter(grepl("dose_", term) & q_value < 0.05) %>%
  dplyr::select(cell_type, term, gene_short_name) %>%
  dplyr::mutate(catalog_number = stringr::str_split_fixed(term, "_", 2)[,2] )

degs_for_plot = rbind(A549_degs_for_plot,K562_degs_for_plot, MCF7_degs_for_plot)

degs_for_plot = left_join(degs_for_plot, drug_annotations, by = "catalog_number")

degs_for_plot %>%
  dplyr::select(-term) %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(x = pathway_level_1, fill = pathway_level_1), size = .2, color = "black") +
  monocle_theme_opts() +
  facet_wrap(~cell_type) +
  xlab("Pathway Targeted") +
  ylab("Number of DEGs") +
  scale_y_log10() +
  coord_flip() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.2,"line"))+
  scale_fill_manual("Pathway", values = pathway_colors) +
ggsave("supplemental_figure_8_DEGs_by_pathway_barplot.pdf", height = 3, width = 6)


