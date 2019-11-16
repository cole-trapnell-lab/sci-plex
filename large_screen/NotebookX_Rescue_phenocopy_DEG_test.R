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
  library(Matrix)
  library(monocle3)
})

# Load processed cds objects from sciCHEM 20
cds.list <- readRDS("A549_MCF7_rescue_phenocopy_processed_cds.rds")

# Perform differential gene expression analysis as a function of progression along the HDACi trajectory
hdac_pseudotime_diff_test_A549 <- fit_models(cds.list[["A549"]], 
                                             model_formula_str = "~splines::ns(Pseudotime, df = 3)",
                                             cores = 40)

hdac_pseudotime_diff_test_A549 <- coefficient_table(hdac_pseudotime_diff_test_A549)

# saveRDS(hdac_pseudotime_diff_test_A549, "hdac_pseudotime_diff_test_A549.rds")

# hdac_pseudotime_diff_test_A549 <- readRDS("hdac_pseudotime_diff_test_A549.rds")

hdac_pseudotime_diff_test_MCF7 <- fit_models(cds.list[["MCF7"]], 
                                             model_formula_str = "~splines::ns(Pseudotime, df = 3)",
                                             cores = 40)

hdac_pseudotime_diff_test_MCF7 <- coefficient_table(hdac_pseudotime_diff_test_MCF7)

# saveRDS(hdac_pseudotime_diff_test_MCF7, "hdac_pseudotime_diff_test_MCF7.rds")

# hdac_pseudotime_diff_test_MCF7 <- readRDS("hdac_pseudotime_diff_test_MCF7.rds")

# Load differential gene expression results for the aligned HDACi trajectory obtained from our large screen

hdac_pseudotime_diff_test_A549 <- readRDS("hdac_pseudotime_diff_test_A549.rds")
hdac_pseudotime_diff_test_MCF7 <- readRDS("hdac_pseudotime_diff_test_MCF7.rds")

hdac_pseudotime_A549_degs <- unique((subset(hdac_pseudotime_diff_test_A549, 
	term %in% c("splines::ns(Pseudotime, df = 3)1","splines::ns(Pseudotime, df = 3)2","splines::ns(Pseudotime, df = 3)3") & 
	q_value < 0.05))$gene_short_name)

hdac_pseudotime_MCF7_degs <- unique((subset(hdac_pseudotime_diff_test_MCF7, 
	term %in% c("splines::ns(Pseudotime, df = 3)1","splines::ns(Pseudotime, df = 3)2","splines::ns(Pseudotime, df = 3)3") & 
	q_value < 0.05))$gene_short_name)

# Load differential gene expression results for the aligned HDACi trajectory obtained from our large screen

hdac_pseudotime_terms <- read.table("hdac_pseudotime_terms.txt", sep = "\t", header = TRUE)
screen_hdac_pseudotime_degs <- unique((subset(hdac_pseudotime_terms, 
	term %in% c("splines::ns(Pseudotime, df = 3)1","splines::ns(Pseudotime, df = 3)2","splines::ns(Pseudotime, df = 3)3") & 
	q_value < 0.05))$gene_short_name)

# Generate a Venn diagram peciting teh overlap between DEG sets

png("Overlap_of_pseudodose_DEGs.png", width = 3, height = 3.25, res = 1000, units = "in")
VennDiagram::draw.triple.venn(length(screen_hdac_pseudotime_degs), 
							  length(hdac_pseudotime_A549_degs), 
							  length(hdac_pseudotime_MCF7_degs), 
							  length(intersect(screen_hdac_pseudotime_degs,hdac_pseudotime_A549_degs)),
							  length(intersect(hdac_pseudotime_A549_degs,hdac_pseudotime_MCF7_degs)), 
							  length(intersect(screen_hdac_pseudotime_degs,hdac_pseudotime_MCF7_degs)),
							  length(unique(Reduce(intersect, 
							  					   list(screen_hdac_pseudotime_degs,
							  					   	    hdac_pseudotime_A549_degs,
							  					   	    hdac_pseudotime_MCF7_degs)))),
							  category = c("Screen","A549","MCF7"),
							  fill = c('red', 'blue', 'green'),
							  euler.d = TRUE,
							  scaled = FALSE,
							  cex = 1,
							  cat.cex = 1,
							  fontfamily = "Helvetica")
dev.off()


