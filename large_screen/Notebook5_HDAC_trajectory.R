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
  library(ggridges)
  library(piano)
  library(tidymodels)
  library(snowfall)
  library(UpSetR)
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



append_umap_coordinates = function(cds){

  num_dim = dim(cds@reducedDims[["UMAP"]])[2]
  for (i in seq(1,num_dim)){
    new_col_name = paste("umap",i,sep = "")
    colData(cds)[[new_col_name]] = cds@reducedDims[["UMAP"]][,i]
  }
  return(cds)
}



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
  filter(!is.na(top_oligo_P),
         !is.na(top_oligo_W)) %>%
  filter(hash_umis_P >= 5,
         hash_umis_W >= 5,
         top_to_second_best_ratio_P >= 5,
         top_to_second_best_ratio_W >= 5) %>%
  filter(rt_well != 1) %>%
  pull(cell)

cds = cds[human.genes,validly_labeled_cells]

# Isolate cells at 24 hours for analysis
cells_24hrs =
  colData(cds) %>%
  as.data.frame() %>%
  filter(time_point == 24) %>%
  pull(cell)

cds = cds[,cells_24hrs]



#################### Do MNN and trajectory reconstruction #####################
hdac_vehicle_cells = read.table("bin/hdac_cells.txt", 
                                header = F, 
                                col.names = "cell")

cds = cds[,colData(cds)$cell %in% hdac_vehicle_cells$cell]


set.seed(seed = 1)

# Down sample cells to equalize representation
cell_sample =
  colData(cds) %>%
  as.data.frame() %>%
  group_by(cell_type) %>%
  # greatest common denominator
  sample_n(size =  12229) %>%
  pull(cell)


cds = cds[,pData(cds)$cell %in% cell_sample]

cds = detect_genes(cds)
cds = estimate_size_factors(cds)
cds <- estimate_cell_cycle(cds,
                           g1s_markers = cc.genes$s.genes,
                           g2m_markers = cc.genes$g2m.genes)

cds <- preprocess_cds(cds,
                      method = 'PCA',
                      num_dim = 25,
                      norm_method = 'log',
                      verbose = T)
#################

cells_list <- list()
cell_types <- unique(colData(cds)[,"cell_type"])

for(genotype in cell_types){
      cells_list[[genotype]] <- row.names(pData(cds)[pData(cds)[,"cell_type"] == genotype,])
}

PCA_cells_1 = cds@reducedDims[["PCA"]][cells_list[[1]],]
PCA_cells_2 = cds@reducedDims[["PCA"]][cells_list[[2]],]
PCA_cells_3 = cds@reducedDims[["PCA"]][cells_list[[3]],]



corrected_PCA <- scran::mnnCorrect(t(PCA_cells_1),
                                        t(PCA_cells_2),
                                        t(PCA_cells_3))

corrected_PCA <- t(do.call("cbind", corrected_PCA[["corrected"]]))

row.names(corrected_PCA) <- c(row.names(PCA_cells_1),
                                   row.names(PCA_cells_2),
                                   row.names(PCA_cells_3))

corrected_PCA <- corrected_PCA[colnames(cds),]
cds@reducedDims[["PCA"]]= as.matrix(corrected_PCA)

################

cds = reduce_dimension(cds,
                       max_components = 2,
                       reduction_method = 'UMAP',
                       umap.metric = 'cosine',
                       umap.n_neighbors = 10,
                       umap.min_dist = 0.1,
                       verbose = TRUE)


cds <- cluster_cells(cds,reduction_method = "PCA")

colData(cds)$Cluster = clusters(cds, reduction_method="PCA")

cds <- cluster_cells(cds,reduction_method = "UMAP")
colData(cds)$louvain_component = cds@clusters[["UMAP"]]$partitions

graph_parameters = list()
graph_parameters[["minimal_branch_len"]] = 10
graph_parameters[["ncenter"]] = 750

cds <- learn_graph(cds, learn_graph_control= graph_parameters, use_partition = T, close_loop = F)
cds = append_umap_coordinates(cds)


# Get the closest vertex for every cell
pData(cds)$closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[,1]

ordering_summary =
  pData(cds) %>%
  as.data.frame() %>%
  dplyr::group_by(closest_vertex) %>%
  dplyr::count(dose) %>%
  dplyr::mutate(total_cells = sum(n), fraction_dose = n / total_cells)

root_nodes =  
  ordering_summary %>% 
  filter(dose == 0 & fraction_dose > 0.5)

root_nodes = root_nodes$closest_vertex

pData(cds)$root_node = pData(cds)$closest_vertex %in% root_nodes
root_cells  = 
  colData(cds) %>% 
  as.data.frame() %>% 
  filter(root_node) %>% pull(cell) %>% 
  as.character()

cds <- order_cells(cds, root_cells=root_cells)

pData(cds)$Pseudotime = cds@principal_graph_aux[["UMAP"]]$pseudotime

pData(cds)$Pseudotime_bin = cut(pData(cds)$Pseudotime, breaks=quantile(pData(cds)$Pseudotime, seq(0, 1, 0.1)), labels=F)

################


colData(cds)$proliferation_cutoff  = 
  ifelse(colData(cds)$cell_type == "MCF7", 2.5, 2)

coldata_cds = 
  colData(cds) %>%
  as.data.frame() %>%
  group_by(cell_type) %>%
  mutate(low_proliferation = proliferation_index < proliferation_cutoff) %>%
  ungroup() %>%
  group_by(Pseudotime_bin,drug_dose) %>%
  add_tally() %>%
  mutate(fraction_in_low = sum(low_proliferation)/n)


coldata_cds %>% 
  filter(!is.na(Pseudotime_bin)) %>%
  dplyr::select(product_name,Pseudotime_bin,fraction_in_low) %>%
  distinct() %>%
  ggplot(aes(x = factor(Pseudotime_bin), y = 100*(fraction_in_low))) +
  geom_boxplot(fill = "grey80", outlier.stroke = 0, outlier.size = 0.5, outlier.colour =  "grey80") +
  geom_jitter(size = 0.5, stroke = 0) +
  monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  xlab("Pseudodose Bin") +
  ylab("Percent Low") +
  ggsave("supplemental_figure29_percent_cells_in_low_proliferation_index.pdf", height = 1.5, width = 2.5)

################

colData(cds)$product_name_alt <- 
  plyr::revalue((pData(cds)$product_name),
                c("Abexinostat (PCI-24781)" = 'Abexinostat', 
                  "Pracinostat (SB939)" = 'Pracinostat',
                  "Mocetinostat (MGCD0103)" = 'Mocetinostat',
                  "Tucidinostat (Chidamide)" = 'Tucidinostat',
                  "AR-42" = 'AR-42', 
                  "Tacedinaline (CI994)" = 'Tacedinaline',
                  "CUDC-101" = 'CUDC-101',
                  "M344" = 'M344',
                  "CUDC-907" = 'CUDC-907',
                  "Dacinostat (LAQ824)" = 'Dacinostat',
                  "Belinostat (PXD101)" = 'Belinostat',
                  "vehicle" = 'Vehicle',
                  "Quisinostat (JNJ-26481585) 2HCl" = 'Quisinostat',
                  "Trichostatin A (TSA)" = 'TSA',
                  "Resminostat" = 'Resminostat',
                  "Entinostat (MS-275)" = 'Entinostat',
                  "Givinostat (ITF2357)" = 'Givinostat',
                  "Panobinostat (LBH589)" = 'Panobinostat',
                  "Vehicle" = 'Vehicle'))

# Figure 4A

colData(cds) %>%
  as.data.frame() %>%
  ggplot()+
  geom_point(aes( x = umap1, y =  umap2, color = dose_character), stroke = 0, size = 0.25)  +
  monocle_theme_opts() +
  scale_color_manual("Log(Dose [nM])",values = c("0" = "gray", "10" = "#1D1147FF",
                                                 "100" = "#822681FF", "1000" = "#E65164FF", "10000" = "#FEC287FF")) +
  theme_void() +
  theme(legend.position = "none", text = element_text(size = 6)) +
  ggsave("figure4_dose_umap_figure4a.png",
         width = 1.5, height = 1.5, unit = "in", dpi = 900)

# Figure 4B,C,D 

pseudotime_bin_summary = 
  pData(cds) %>% 
  as.data.frame() %>%
  group_by(Pseudotime_bin) %>%
  summarise(num_in_bin = n())

pseudotime_bin_summary =
  left_join(as.data.frame(pData(cds)),
            pseudotime_bin_summary, 
            by = "Pseudotime_bin")

pseudotime_bin_summary = 
  pseudotime_bin_summary %>%
  group_by(Pseudotime_bin, cell_type) %>%
  add_tally() %>%
  mutate(fraction_cell_type = n/num_in_bin) %>%
  dplyr::select(-n) %>%
  ungroup()

pseudotime_bin_summary = 
  pseudotime_bin_summary %>%
  group_by(Pseudotime_bin, product_name_alt) %>%
  add_tally() %>%
  mutate(fraction_product_name = n/num_in_bin) %>%
  dplyr::select(-n) %>%
  ungroup()

pseudotime_bin_summary = 
  pseudotime_bin_summary %>%
  group_by(Pseudotime_bin, dose_character) %>%
  add_tally() %>%
  mutate(fraction_dose = n/num_in_bin) %>%
  dplyr::select(-n) %>%
  ungroup()



prop_prod_fill_settings =   
  scale_fill_manual("Product Name", values = c("Panobinostat"="aquamarine3", 
                                               "Resminostat"="lemonchiffon4", 
                                               "Givinostat"="deepskyblue2", 
                                               "Quisinostat"="slategray4", 
                                               "Dacinostat"="navy", 
                                               "CUDC-907"="brown4", 
                                               "M344"="darkgreen", 
                                               "CUDC-101"="orangered2", 
                                               "Tacedinaline" = "gold2", 
                                               "AR-42" = "darkolivegreen",
                                               "Tucidinostat" = "chartreuse3",
                                               "Pracinostat" = "plum4", 
                                               "Mocetinostat" = "darkorchid2",
                                               "Abexinostat" = "darkcyan",
                                               "Entinostat" = "firebrick1",
                                               "Belinostat" ="darkgoldenrod4",
                                               "TSA" = "springgreen1",
                                               "Vehicle" = "grey80"))


pseudotime_bin_summary %>%
  dplyr::select(Pseudotime_bin,dose_character, fraction_dose ) %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(x = Pseudotime_bin, y =fraction_dose, fill = dose_character ), 
           color = "black", size = .25, stat = "identity") +
  scale_fill_manual("Dose (nM)",
                    values = c("0" = "gray", "10" = "#1D1147FF","100" = "#822681FF", 
                               "1000" = "#E65164FF", "10000" = "#FEC287FF")) + 
  monocle:::monocle_theme_opts() +
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +    
  scale_x_continuous(breaks = seq(1,10,1))+
  xlab("Pseudodose Bin") +
  ylab("Proportion") +
  ggsave("supplemental_figure21_dose_proportions.pdf",
         width = 2, height = 1.5, unit = "in", dpi = 1200)

pseudotime_bin_summary %>%
  dplyr::select(Pseudotime_bin,product_name_alt, fraction_product_name ) %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(x = Pseudotime_bin, y =fraction_product_name, fill = product_name_alt ), 
           color = "black", size = .25, stat = "identity") +
  monocle:::monocle_theme_opts() +
  prop_prod_fill_settings +
  scale_x_continuous(breaks = seq(1,10,1))+
  guides(fill=guide_legend(ncol=2)) + 
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +  
  xlab("Pseudodose Bin") +
  ylab("Proportion") +
  ggsave("supplemental_figure21_product_proportions.pdf",
         width = 2.5, height = 1.5, unit = "in", dpi = 1200)


#########
ica_space_df <- 
  t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>%
  as.data.frame() %>% dplyr::select(prin_graph_dim_1 = V1,
                                     prin_graph_dim_2 = V2) %>% 
  dplyr::mutate(sample_name = rownames(.),sample_state = rownames(.))

edge_df <- 
  cds@principal_graph[["UMAP"]] %>% 
  igraph::as_data_frame() %>% 
  dplyr::select_(source = "from",target = "to") %>% 
  dplyr::left_join(ica_space_df %>% 
                     dplyr::select(source = sample_name,
                                   source_prin_graph_dim_1 = prin_graph_dim_1,
                                   source_prin_graph_dim_2 = prin_graph_dim_2),
                   by = "source") %>% 
  dplyr::left_join(ica_space_df %>%
                   dplyr::select(target = sample_name, 
                                   target_prin_graph_dim_1 = prin_graph_dim_1,
                                   target_prin_graph_dim_2 = prin_graph_dim_2),
                   by = "target")


ggplot() +
  geom_point(data = 
               colData(cds) %>%
               as.data.frame(),
             aes(x = -umap1, y = -umap2,color = dose_character), 
             size = .25, stroke = 0)+
  geom_segment(data = edge_df,
               aes_string(x="source_prin_graph_dim_1", 
                          y="source_prin_graph_dim_2", 
                          xend="target_prin_graph_dim_1", 
                          yend="target_prin_graph_dim_2"), 
               size=.75, linetype="solid", na.rm=TRUE) +
  geom_point(data = colData(cds) %>% 
               as.data.frame() %>%
               filter(root_node) %>%
               dplyr::select(source = closest_vertex) %>%
               mutate(source = paste("Y_", source, sep = "")) %>%
               left_join(edge_df, by = "source"), 
             aes(x = source_prin_graph_dim_1, y = source_prin_graph_dim_2), 
             color = "red", size = 1.5, stroke = 0)+ 
  scale_color_manual(values = c("0" = "gray", "10" = "#1D1147FF",
                                "100" = "#822681FF", "1000" = "#E65164FF", 
                                "10000" = "#FEC287FF")) +
  theme_void()+
  theme(legend.position = "none", text = element_text(size = 6)) +
  ggsave("supplemental_figure21_UMAP_root_nodes.png",
         width = 1.5, height = 1.5, unit = "in", dpi = 1200)  

pseudotime_bin_summary %>%
  dplyr::select(Pseudotime_bin,cell_type, fraction_cell_type ) %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(x = Pseudotime_bin, y =fraction_cell_type, fill = cell_type ), 
           color = "black", size = .15, stat = "identity") +
  scale_fill_manual("Cell Line",values = c("A549" = "#438FCD","K562" = "#D6C939","MCF7" ="#DB3C6A")) +
  monocle:::monocle_theme_opts() +
  scale_x_continuous(breaks = seq(1,10,1))+
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +    
  xlab("Pseudodose Bin") +
  ylab("Proportion") +
  ggsave("supplemental_figure21_cell_type_proportions_bins.pdf",
         width = 2, height = 1.5, unit = "in", dpi = 1200)

colData(cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(data = 
               pData(cds) %>%
               as.data.frame() %>%
               dplyr::select(-cell_type),
             aes(x = umap1, y = umap2), 
             stroke = 0, size =.25, color = "grey80", alpha = .5) +
  geom_density_2d(aes(x = umap1, y = umap2, color = cell_type),size =.25) +
  theme_void()+
  scale_color_manual(values = c("#438FCD","#D6C939","#DB3C6A" )) +
  facet_wrap(~cell_type) +
  theme(legend.position = "none", 
        text = element_text(size = 6),
        strip.text.x = element_text(size = 6)) + 
  ggsave("supplemental_figure21_cell_type_densities.png",
         height = 1.5, width = 3, dpi = 1200)

vehicle.df = 
  pData(cds) %>% 
  as.data.frame() %>%
  dplyr::select(Pseudotime, dose_character, product_name_alt, cell_type) %>% 
  filter(product_name_alt == "Vehicle") %>% 
  dplyr::select(-product_name_alt)

non.vehicle.df = 
  pData(cds)  %>% 
  as.data.frame() %>%
  dplyr::select(Pseudotime, dose_character, product_name_alt, cell_type) %>% 
  filter(product_name_alt != "Vehicle")

unique.keys.df = 
  pData(cds)  %>% 
  as.data.frame() %>%
  dplyr::select(cell_type, product_name_alt) %>% 
  unique()

duplicated.vehicle.df = inner_join(unique.keys.df, vehicle.df, by = "cell_type")

duplicated.vehicle.df = 
  duplicated.vehicle.df %>%
  filter(product_name_alt != "Vehicle") %>%
  dplyr::select(cell_type, product_name_alt, Pseudotime, dose_character)

non.vehicle.df = 
  non.vehicle.df %>%
  dplyr::select(cell_type, product_name_alt, Pseudotime, dose_character)

combined.df = rbind(duplicated.vehicle.df,non.vehicle.df)

combined.df %>%
  filter(is.finite(Pseudotime)) %>%
  ggplot() +
  geom_density_ridges(aes( x = Pseudotime, y = dose_character, fill = dose_character), size = .25) + 
  monocle:::monocle_theme_opts() +
  facet_wrap(product_name_alt~cell_type,ncol = 6) +
  scale_fill_manual("Log(Dose [nM])",values = c("0" = "gray", "10" = "#1D1147FF",
                                                "100" = "#822681FF", "1000" = "#E65164FF", "10000" = "#FEC287FF")) +
  theme(legend.position = "none", 
        text = element_text(size = 6),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("") +
  ylab("")+
  ggsave("supplemental_figure22_ridges_hdac_trajectory.pdf", height = 8 , width = 6, unit = "in")



combined.df %>%
  filter(is.finite(Pseudotime)) %>%
  filter(product_name_alt == "Pracinostat") %>%
  ggplot() +
  geom_density_ridges(aes( x = Pseudotime, y = dose_character, fill = dose_character), size = .15) + 
  monocle:::monocle_theme_opts() +
  facet_wrap(~cell_type,ncol = 3) +
  scale_fill_manual("Log(Dose [nM])",values = c("0" = "gray", "10" = "#1D1147FF",
                                                "100" = "#822681FF", "1000" = "#E65164FF", "10000" = "#FEC287FF")) +
  theme(legend.position = "none", 
        text = element_text(size = 6),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_blank()) + 
  xlab("") +
  ylab("")+
  ggsave("figure4_Fig4_ridges_Pracinostat.pdf", height = .55 , width = 1.5, unit = "in")

combined.df %>%
  filter(is.finite(Pseudotime)) %>%
  filter(product_name_alt == "Quisinostat") %>%
  ggplot() +
  geom_density_ridges(aes( x = Pseudotime, y = dose_character, fill = dose_character), size = .15) + 
  monocle:::monocle_theme_opts() +
  facet_wrap(~cell_type,ncol = 3) +
  scale_fill_manual("Log(Dose [nM])",values = c("0" = "gray", "10" = "#1D1147FF",
                                                "100" = "#822681FF", "1000" = "#E65164FF", "10000" = "#FEC287FF")) +
  theme(legend.position = "none", 
        text = element_text(size = 6),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_blank()) +   
  xlab("") +
  ylab("")+
  ggsave("figure4_Fig4_ridges_Quis.pdf", height = .55 , width = 1.5, unit = "in")

combined.df %>%
  filter(is.finite(Pseudotime)) %>%
  filter(product_name_alt == "CUDC-907") %>%
  ggplot() +
  geom_density_ridges(aes( x = Pseudotime, y = dose_character, fill = dose_character), size = .15) + 
  monocle:::monocle_theme_opts() +
  facet_wrap(~cell_type,ncol = 3) +
  scale_fill_manual("Log(Dose [nM])",values = c("0" = "gray", "10" = "#1D1147FF",
                                                "100" = "#822681FF", "1000" = "#E65164FF", "10000" = "#FEC287FF")) +
  theme(legend.position = "none", 
        text = element_text(size = 6),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_blank()) +   
  xlab("") +
  ylab("")+
  ggsave("figure4_ridges_CUDC907.pdf", height = .55 , width = 1.5, unit = "in")


saveRDS(object = cds, file = "hdac_aligned_cds.RDS")
