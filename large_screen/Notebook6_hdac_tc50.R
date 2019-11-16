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
  library(stringr)
  library(tidyr)
  library(ggridges)
  library(devtools)
  library(tidymodels)
  library(drc)
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
                "dose_response_TC50.R",
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

cds_coldata = 
  colData(cds) %>%
  as.data.frame()

cds_coldata$treatment  = cds_coldata$product_name_alt

sci_Plex_pseudodoses = measure_pseudodose(cds_coldata)

sci_Plex_models = fit_dose_response_models(sci_Plex_pseudodoses, max_pseudodose=max(sci_Plex_pseudodoses$pseudodose), verbose=T)

sci_Plex_model_stats = extract_dose_model_stats(sci_Plex_pseudodoses, sci_Plex_models)

write.table(sci_Plex_model_stats, "Supplementary_Table_6_TC50_models.txt",sep = "\t",quote = F )

pdf("supplemental_figure27_HDAC_TC50_curves_all.pdf", width=6, height=24)
plot_dose_response_curves(sci_Plex_pseudodoses, 
                          sci_Plex_models, ncol=3) + theme(text = element_text(size = 6)) 
dev.off()

drugs_to_plot = c("Belinostat", "TSA")

pdf("supplemental_figure27_select_TC50_curves.pdf", width=4, height=2)
plot_dose_response_curves(sci_Plex_pseudodoses %>% filter(treatment %in% drugs_to_plot), 
                          sci_Plex_models %>% filter(treatment %in% drugs_to_plot), ncol=3) + 
  theme(text = element_text(size = 6)) 
dev.off()


sci_Plex_model_stats_anot = 
  sci_Plex_model_stats %>% 
  dplyr::select(everything(), product_name_alt = treatment)

prod_names = cds_coldata %>% dplyr::select(catalog_number, treatment, product_name_alt) %>% distinct()

sci_Plex_model_stats_anot = 
  left_join(sci_Plex_model_stats_anot,prod_names, by = "product_name_alt" )


ic50_table = tibble::as_tibble(read.table(file = "bin/hdacic50.txt", header= T, sep = ","))
ic50_matrix = as.matrix(ic50_table[,-c(1,2)])
ic50_matrix[ic50_matrix >= 10000] = NA
ic50_matrix = ic50_matrix/1000000000
aggregegate_ic50 = rowMeans(log10(ic50_matrix), na.rm=TRUE)
ic50_table$aggregate_ic50 = aggregegate_ic50
ic50_table = inner_join(ic50_table, sci_Plex_model_stats_anot, by="catalog_number")

solubility_table = read.table("bin/compound_solubility.csv", header = F, 
                              col.names = c("catalog_number","solubility"), sep = ",")

ic50_table = left_join(ic50_table, solubility_table, by="catalog_number")


pdf("figure4_HDAC_TC50_vs_IC50.pdf", width=3, height=2.25)
ggplot(ic50_table) +
  geom_smooth(data = ic50_table %>% 
                filter(!(product_name_short %in% c("Mocetinostat", "Trichostatin", "CUDC-101"))), 
              aes(x = aggregate_ic50 ,y = log(ec50)), se=F) + 
  geom_point(data = ic50_table %>% filter(!(product_name_short %in% c("Mocetinostat", "Trichostatin", "CUDC-101"))), 
             aes(x = aggregate_ic50 ,y = log(ec50), color = product_name_alt), 
             stroke = 0, size =1.75)+
  geom_point(data = ic50_table %>% filter((product_name_short %in% c("Mocetinostat", "Trichostatin", "CUDC-101"))), 
             aes(x = aggregate_ic50 ,y = log(ec50), color = product_name_alt), 
             size =5, shape = 42, show.legend = F)+
  
  scale_color_manual("Product Name", values = c("Panobinostat"="aquamarine3", 
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
                                                "Vehicle" = "grey80"))+
  
  theme(legend.position = "right", 
        text = element_text(size = 7),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line"))+      
  monocle3:::monocle_theme_opts() +
  xlab("mean(log10(IC50))") +
  ylab("log(TC50)")

dev.off()


pdf("supplemental_figure27_HDAC_TC50_vs_IC50_solubility.pdf", width=2.5, height=2.1)
ggplot(ic50_table) +
  geom_point(data = ic50_table %>% filter(!(product_name_short %in% c("Mocetinostat", "Trichostatin", "CUDC-101"))),
             aes(x = aggregate_ic50 ,y = log(ec50), color = solubility), 
             stroke = 0, size =1.75)+
  geom_point(data = ic50_table %>% filter((product_name_short %in% c("Mocetinostat", "Trichostatin", "CUDC-101"))), 
             aes(x = aggregate_ic50 ,y = log(ec50), color = solubility), 
             size =4, shape = 42, show.legend = F)+
  scale_color_viridis_c() +
  guides(color = guide_colourbar(barwidth = 0.5, barheight = 3)) +
  theme(legend.position = "right", 
        text = element_text(size = 7))+      
  monocle3:::monocle_theme_opts() +
  xlab("mean(log10(IC50))") +
  ylab("log(TC50)")

dev.off()


ic50_by_hdac_iso = 
  ic50_table %>%
  dplyr::select(starts_with("ic50"),cell_type,product_name_alt,solubility,ec50 ) %>% 
  distinct() %>%
  gather(starts_with("IC50_H"), key = enzyme, value = ic50) %>%
  mutate(ic50 = ic50/1000000000) %>%
  separate(col = enzyme, into = c("dummy","HDAC"), sep = "_") %>%
  filter(!(product_name_alt %in% c("Mocetinostat", "TSA", "CUDC-101")))


pdf("supplemental_figure27_hdac_potentcy_facet.pdf", width=6, height=2.75)
ggplot(ic50_by_hdac_iso, aes(x = log10(ic50) ,y = log(ec50))) +
  geom_point(aes(color = product_name_alt), 
             stroke = 0, size =1.75)+
  monocle:::monocle_theme_opts() +
  facet_wrap(~HDAC,scales = "free", ncol = 4) +
  theme(legend.position = "right", 
        text = element_text(size = 7),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +  
  scale_color_manual("Product Name", values = c("Panobinostat"="aquamarine3", 
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
                                                "Vehicle" = "grey80")) +
  xlab("log10(IC50)") +
  ylab("log(TC50)")
dev.off()
