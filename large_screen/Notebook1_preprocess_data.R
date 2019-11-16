# devtools::install_github(repo = 'cole-trapnell-lab/monocle3', ref = "d4a9a35")
# Specific commit used in performing all analyses

# Input path to github directory
path_to_github = "~/sci-Plex/"

# Set directory for sciPlex github 
bin_directory = paste0(path_to_github,
                       "bin",
                       sep = "")

suppressPackageStartupMessages({
  library(dplyr)
  library(monocle3)
})


DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)

setwd(paste0(path_to_github,
             "large_screen",
             sep = ""))


# Read in the CDS
cds = readRDS("sciPlex3_cds.RDS")


################## Append Cell Path Metadata To CDS ################## 
#Pull out CDS column data to append relevant metadata
coldata_cds = colData(cds) %>% as.data.frame()

coldata_cds = 
  coldata_cds %>%
  dplyr::rename(cell = Cell,
                rt_well = RTwell,
                lig_well = Ligwell,
                pcr_well = PCR_well,
                pcr_plate = PCR_plate) %>%
  dplyr::mutate(culture_plate = top_oligo_P)


rt_plate = 
  data.frame(rt_well = seq(1,384),
             rt_plate = seq(0,383) %/% 96  + 1)
lig_plate = 
  data.frame(lig_well = seq(1,384),
             lig_plate = seq(0,383) %/% 96  + 1)

coldata_cds = left_join(coldata_cds, rt_plate, by = "rt_well")
coldata_cds = left_join(coldata_cds, lig_plate, by = "lig_well")

coldata_cds$Combo = paste0(coldata_cds$top_oligo_W,coldata_cds$top_oligo_P)

################## Append Hash Encoding Metadata To CDS ################## 

hash_meta_data = read.table("bin/sciChem3_screen2_hashMeta_data.tsv",
                            header = F, 
                            sep = "\t")

colnames(hash_meta_data) = 
  c("well_oligo","plate_oligo","cell_type","replicate","time_point","drug_dose")

hash_meta_data$Combo =  paste0(hash_meta_data$well_oligo,hash_meta_data$plate_oligo)
coldata_cds = left_join(coldata_cds, hash_meta_data, by = "Combo")

coldata_cds$catalog_number = stringr::str_split_fixed(coldata_cds$drug_dose, pattern = "_",n = 2)[,1]

coldata_cds$vehicle = grepl(pattern = "S000",x = coldata_cds$catalog_number)

coldata_cds$dose_pattern = stringr::str_split_fixed(coldata_cds$drug_dose, pattern = "_",n = 2)[,2]

dose_keys = data.frame(dose_pattern = c("1","2","3","4"),
           dose_character = c("10000","1000","100","10"))

coldata_cds = left_join(coldata_cds,dose_keys, by = "dose_pattern")

# Set the dose of vehicle cells to 0
coldata_cds$dose_character = ifelse(coldata_cds$vehicle, "0", as.character(coldata_cds$dose_character))
coldata_cds$dose = as.numeric(coldata_cds$dose_character)

coldata_cds$catalog_number = ifelse(coldata_cds$vehicle, "S0000", as.character(coldata_cds$catalog_number))
coldata_cds$treatment = as.character(coldata_cds$catalog_number)

################## Append Drug Pathway Metadata To CDS ################## 

pathway_annotations = read.table("bin/Supplementary_Table_3.txt",
                                 header = T, sep = "\t")

pathway_annotations = 
  pathway_annotations %>%
  dplyr::select(pathway_level_1,
                pathway_level_2,
                catalog_number) %>%
  dplyr::mutate(catalog_number  = as.character(catalog_number)) %>%
  dplyr::distinct()


pathway_annotations$catalog_number[is.na(pathway_annotations$catalog_number)] <- "S0000"

coldata_cds = left_join(coldata_cds,pathway_annotations, by ="catalog_number" )

################## Append Drug Metadata To CDS ################## 

drug_annotations = read.table("bin/drugProperties_short.tsv",
                              header = T, 
                              sep = "\t")

drug_annotations =
  drug_annotations %>%
  dplyr::rename(catalog_number = Catalog.Number,
                product_name = Product.Name,
                pathway = Pathway,
                target = Target)


coldata_cds = left_join(coldata_cds,drug_annotations, by = "catalog_number")


################## Write and save new CDS ################## 

rownames(coldata_cds) = coldata_cds$cell
mat = counts(cds)
gene_metadata = rowData(cds)

cds = new_cell_data_set(mat,
                        cell_metadata =coldata_cds,
                        gene_metadata = gene_metadata)

saveRDS(object = cds, 
        file = "cds.RDS")

