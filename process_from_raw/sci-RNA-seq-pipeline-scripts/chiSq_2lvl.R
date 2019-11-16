suppressPackageStartupMessages({
  library(ggplot2)
  library(dirmult)
  library(tidyr)
  library(dplyr)
  library(parallel)
  library(pbapply)
  library(UpSetR)
  library(monocle)
})


args = commandArgs(trailingOnly = T)

if (length(args) < 5) {
    cat("Usage: Rscript chiSqTest.R background_cells_list int_cells_list real_cells_list hashTable_Path", file = stderr())
    exit(1)
}

background_cells_list = args[1]
real_cells_list = args[2]
hashed_fastq_path = as.character(args[3])
hashSample_sheet_path = args[4]
cds_path = args[5]
output.dir = args[6]

prepare_hash_matrix <- function(hashes){
    hash_matrix = 
        hashes %>% 
        select(Cell, Oligo, Count) %>% 
        distinct() %>% 
        spread(key="Oligo", value="Count", fill = 0, drop = F)
    row.names(hash_matrix) = hash_matrix$Cell
    cell_ids = hash_matrix$Cell 
    hash_matrix = hash_matrix[,-1]
    hash_matrix = as.matrix(hash_matrix)

    return(list(hash_matrix,cell_ids))
}

generate_pval_matrix <- function(test_hash_matrix, hash_frequencies, overdispersion, cores=1){
    # Note: use outfile argument to makeCluster for debugging
    platform <- Sys.info()[['sysname']]
    if (platform == "Windows")
        cl <- makeCluster(cores)
    if (platform %in% c("Linux", "Darwin")) 
        cl <- makeCluster(cores)
    
    cleanup <- function(){
        stopCluster(cl)
    }
    on.exit(cleanup)
    required_packages = c("dirmult")
    if (is.null(required_packages) == FALSE){
        clusterCall(cl, function(pkgs) {
        for (req in pkgs) {
            library(req, character.only=TRUE)
        }
        }, required_packages)
    }
    pval_matrix = t(pbapply(test_hash_matrix, 1, function(x, hash_frequencies, overdispersion) {
        total_counts = sum(x)
        sim_sample = simPop(J=1000, 
            n=total_counts,
            pi=hash_frequencies, 
            theta=overdispersion)
        cdf_list = apply(sim_sample$data, 2, function(y) { ecdf(y) })
        pvals = list()
        for (i in 1:length(cdf_list)){
            pvals = append(pvals, cdf_list[[i]] (x[i]) )
        }
        pvals = 1 - unlist(pvals)
        qvals = p.adjust(pvals)
        qvals
    }, 
    hash_frequencies=background_hash_frequencies, 
    overdispersion=overdispersion, 
    cl=cl))
    colnames(pval_matrix) = colnames(test_hash_matrix)
    return(pval_matrix)
}


fit_background_dist <- function(hash_matrix, training_samples=5000, return_training_matrix=FALSE){
    #hash_mtrix = hash_matrix[sample(nrow(hash_matrix), 100),]
    hash_mtrix_sample = hash_matrix[sample(nrow(hash_matrix), training_samples),]
    background_fit = dirmult(hash_mtrix_sample, trace=TRUE)
    #dirmult.summary(hash_mtrix_sample, background_fit)
    if (return_training_matrix)
        return(list(background_fit, hash_mtrix_sample))
    else
        return(background_fit)
}


chisq_vs_background <- function(test_hash_matrix, hash_frequencies){
    hash_frequencies_nz = which(hash_frequencies > 0)
    hash_frequencies = hash_frequencies[hash_frequencies_nz]
    pvals= pbapply(test_hash_matrix[,hash_frequencies_nz], 1, function(x) {
        tryCatch({
            res = chisq.test(x, p=hash_frequencies,  simulate.p.value = FALSE)
            unlist(res[["p.value"]])
            }, error = function(e) { 1.0 }
        )
    })
    return(pvals)
}



assign_hash_labels <- function(test_cell_hashes, background_cell_hashes, min_best_vs_second_best_ratio=2, qval_thresh = 0.05, downsample_rate=NULL){
  background_hash_matrix = prepare_hash_matrix(background_cell_hashes)[[1]]
  # if (is.null(downsample_rate) == FALSE){
  #     background_hash_matrix = floor(background_hash_matrix * downsample_rate)
  # }
  background_hash_matrix = background_hash_matrix[rowSums(background_hash_matrix)>0,]
  
  background_hash_frequencies = colSums(background_hash_matrix)/sum(colSums(background_hash_matrix))
  
  test_hash_list = prepare_hash_matrix(test_cell_hashes)
  test_hash_matrix = test_hash_list[[1]]
  if (is.null(downsample_rate) == FALSE){
    test_hash_matrix = floor(test_hash_matrix * downsample_rate)
  }
  pvals = chisq_vs_background(test_hash_matrix, hash_frequencies=background_hash_frequencies)
  qvals = p.adjust(pvals)                             
  
  expected_background_hashes = outer(rowSums(test_hash_matrix), background_hash_frequencies)
  background_subtracted_test_hashes = test_hash_matrix - expected_background_hashes
  background_subtracted_test_hashes[background_subtracted_test_hashes < 0] = 0
  #hash_hits = background_subtracted_test_hashes * (qvals < qval_thresh)
  top_to_second_best_ratios = apply(background_subtracted_test_hashes, 1, function(x) { y = sort(x, decreasing=TRUE); y[1] / y[2]})
  top_hash = apply(background_subtracted_test_hashes, 1, function(x) { 
    m = which(x == max(x));
    if (length(m) > 1)
      return (NA)
    colnames(background_subtracted_test_hashes)[which(x == max(x))] 
  })
  #unambiguous_hits = which(top_to_second_best_ratios > min_best_vs_second_best_ratio)
  hash_label_df = data.frame(Cell = test_hash_list[[2]],
                             hash_umis = rowSums(test_hash_matrix),
                             pval = pvals,
                             qval = qvals,
                             #hit = qvals < qval_thresh, 
                             top_to_second_best_ratio = top_to_second_best_ratios,
                             #unambiguous_hit = top_to_second_best_ratios > min_best_vs_second_best_ratio,
                             top_oligo = top_hash)
}

#########

# Read in hash Oligos

hashSample_sheet = read.table(hashSample_sheet_path, header = F)
colnames(hashSample_sheet) = c("Oligo", "Seq", "Axis")

# Read in gzipped hash files
hash.files = list.files(path = hashed_fastq_path, pattern ="*.gz")
# Iterate through each file while removing duplicates 
# and counting the number of UMIs per hash
hash_dfs  =lapply(hash.files, function(f) {
  print(paste0("Reading in ",hashed_fastq_path,f))
  tmp.df = read.delim(gzfile(paste0(hashed_fastq_path,f)), header = F) 
  colnames(tmp.df) = c("SampleName", "Cell", "UMI", "Oligo", "Axis")
  tmp.df$Cell = as.character(tmp.df$Cell)
  tmp.df$Oligo = factor(tmp.df$Oligo,levels =  hashSample_sheet$Oligo)
  tmp.df %>% 
    distinct() %>%
    group_by(SampleName,Cell,Axis,Oligo) %>%
    summarise(Count = n()) %>%
    ungroup()
  
})

# Read in list of background cells 
background_cells = read.csv(background_cells_list)
background_cells = as.character(background_cells[,1])


# Read in list of real cells 
real_test_cell_df = read.csv(real_cells_list)
test_cells = as.character(real_test_cell_df[,1])


# Set up containers for result 
hash_df_all <- data.frame(hash_umis=as.numeric(),
                               pval=numeric(), 
                               qval=numeric(),
                               top_to_second_best_ratio = numeric(),
                               top_oligo=character()) 



downsample_rate =1
#background_sample_size = 5000

for( i in 1:length(hash_dfs) ) {
   
  print(paste0("Starting hash assignment ", i))

  # Get the cells that are in the backgroud for this PCR well
  background_cells_curr = background_cells[background_cells %in% hash_dfs[[i]]$Cell]
  # Sample 10K cells to ascertain the background distribution
  #background_cell_sample = background_cells_curr[sample(length(background_cells_curr), background_sample_size)]
  
  background_sample_hashes = hash_dfs[[i]] %>% filter(Cell %in% background_cells_curr, Axis == 1)
  
  
  test_cell_hashes = subset(hash_dfs[[i]], Cell %in% test_cells & Axis == 1)
  hash_df_tmp = assign_hash_labels(test_cell_hashes, background_sample_hashes, downsample_rate=downsample_rate)
  hash_df_all = rbind(hash_df_all, hash_df_tmp)
  
}

# Clean up the result 

cds = readRDS(cds_path)
pData(cds)$Cell = as.character(pData(cds)$Cell)

pData(cds) = left_join(pData(cds), hash_df_all,  by = "Cell")
rownames(pData(cds)) = colnames(exprs(cds))




saveRDS(object = cds,
        file = paste(output.dir, "/", "cds.RDS", sep = ""))




