
collect_gsa_hyper_results_clusters <- function (genes_list, clusters, gsc) 
{
    gene_universe <- unique(as.character(genes_list))
    gsa_results <- list()
    cluster_ids <- unique(clusters)
    for (i in (1:length(cluster_ids))) {
        cluster_genes <- unique(names(clusters[clusters == i]))
        gsaRes <- runGSAhyper(cluster_genes, gsc = gsc, universe = gene_universe, 
            adjMethod = "BH")
        gsa_results[[length(gsa_results) + 1]] <- gsaRes
    }
    names(gsa_results) <- cluster_ids
    gsa_results
}
                                   
gsea_bar_plots <- function(GSAhyper_list, qval_cutoff, pattern, width, height, sample, gsc){
    
    for(cluster in names(GSAhyper_list)){
        
        print(cluster)
        
            GSAhyper_df <- as.data.frame(GSAhyper_list[[cluster]]$p.adj)
            GSAhyper_df$gene_set <- row.names(GSAhyper_df)
            colnames(GSAhyper_df) <- c("qval","gene_set")

            if(is.null(pattern) == FALSE){
                GSAhyper_df$gene_set <- stringr::str_replace(string = GSAhyper_df$gene_set, pattern = pattern, replace = "")
            }

            GSAhyper_df_cutoff <- GSAhyper_df %>% filter(qval < qval_cutoff) %>% arrange(desc(qval)) %>% 
            mutate(gene_set = factor(gene_set, levels = gene_set))

            plot_title <- paste0(sample,"_",as.character(cluster),"_",gsc,".png")
            print(plot_title)
            
            ggplot(GSAhyper_df_cutoff, aes(x = gene_set, y = -log10(qval))) + 
            geom_bar(stat = "identity") + 
            coord_flip() +
            theme_classic(base_size = 8) +
            ggsave(plot_title, width = width, height = height)
        
        
        
    }
    
}

loadGSCSafe <- function (file, type = "auto", addInfo, sep="\t", encoding="latin1") 
{
  if (missing(addInfo)) {
    addUserInfo <- "skip"
    addInfo <- "none"
  }
  else {
    addUserInfo <- "yes"
  }
  tmp <- try(type <- match.arg(type, c("auto", "gmt", "sbml", 
                                       "sif", "data.frame"), several.ok = FALSE), silent = TRUE)
  if (class(tmp) == "try-error") {
    stop("argument type set to unknown value")
  }
  if (type == "auto") {
    if (class(file) == "character") {
      tmp <- unlist(strsplit(file, "\\."))
      type <- tolower(tmp[length(tmp)])
      if (!type %in% c("gmt", "sif", "sbml", "xml")) 
        stop(paste("can not handle .", type, " file extension, read manually using e.g. read.delim() and load as data.frame", 
                   sep = ""))
    }
    else {
      type <- "data.frame"
    }
  }
  if (type == "gmt") {
    con <- file(file, encoding=encoding)
    tmp <- try(suppressWarnings(open(con)), silent = TRUE)
    if (class(tmp) == "try-error") 
      stop("file could not be read")
    if (addUserInfo == "skip") 
      addInfo <- vector()
    gscList <- list()
    i <- 1
    tmp <- try(suppressWarnings(while (length(l <- scan(con, 
                                                        nlines = 1, what = "character", quiet = T, sep=sep)) > 0) {
      if (addUserInfo == "skip") 
        addInfo <- rbind(addInfo, l[1:2])
      tmp <- l[3:length(l)]
      gscList[[l[1]]] <- unique(tmp[tmp != "" & tmp != 
                                      " " & !is.na(tmp)])
      i <- i + 1
    }), silent = TRUE)
    if (class(tmp) == "try-error") 
      stop("file could not be read")
    close(con)
    gsc <- gscList[!duplicated(names(gscList))]
    if (addUserInfo == "skip") 
      addInfo <- unique(addInfo)
  }
  else if (type %in% c("sbml", "xml")) {
    require(rsbml)
    tmp <- try(sbml <- rsbml_read(file))
    if (class(tmp) == "try-error") {
      stop("file could not be read by rsbml_read()")
    }
    gsc <- list()
    for (iReaction in 1:length(reactions(model(sbml)))) {
      metIDs <- names(c(reactants(reactions(model(sbml))[[iReaction]]), 
                        products(reactions(model(sbml))[[iReaction]])))
      geneIDs <- names(modifiers(reactions(model(sbml))[[iReaction]]))
      if (length(geneIDs) > 0) {
        geneNames <- rep(NA, length(geneIDs))
        for (iGene in 1:length(geneIDs)) {
          geneNames[iGene] <- name(species(model(sbml))[[geneIDs[iGene]]])
        }
        for (iMet in 1:length(metIDs)) {
          gsc[[metIDs[iMet]]] <- c(gsc[[metIDs[iMet]]], 
                                   geneNames)
        }
      }
    }
    if (length(gsc) == 0) {
      stop("no gene association found")
    }
    else {
      for (iMet in 1:length(gsc)) {
        tmp1 <- name(species(model(sbml))[[names(gsc)[iMet]]])
        tmp2 <- compartment(species(model(sbml))[[names(gsc)[iMet]]])
        names(gsc)[iMet] <- paste(tmp1, " (", tmp2, ")", 
                                  sep = "")
      }
    }
  }
  else if (type == "sif") {
    tmp <- try(gsc <- as.data.frame(read.delim(file, header = FALSE, 
                                               quote = "", as.is = TRUE), stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("argument file could not be read and converted into a data.frame")
    }
    if (ncol(gsc) != 3) {
      stop("sif file should contain three columns")
    }
    if (addUserInfo == "skip") 
      addInfo <- gsc[, c(1, 2)]
    gsc <- gsc[, c(3, 1)]
    tmp <- nrow(gsc)
    gsc <- unique(gsc)
    geneSets <- unique(gsc[, 2])
    gscList <- list()
    for (iGeneSet in 1:length(geneSets)) {
      gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet], 
                                 1]
    }
    names(gscList) <- geneSets
    gsc <- gscList
  }
  else if (type == "data.frame") {
    tmp <- try(gsc <- as.data.frame(file, stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("argument file could not be converted into a data.frame")
    }
    for (i in 1:ncol(gsc)) {
      gsc[, i] <- as.character(gsc[, i])
    }
    if (ncol(gsc) != 2) {
      stop("argument file has to contain exactly two columns")
    }
    tmp <- nrow(gsc)
    gsc <- unique(gsc)
    geneSets <- unique(gsc[, 2])
    gscList <- list()
    for (iGeneSet in 1:length(geneSets)) {
      gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet], 
                                 1]
    }
    names(gscList) <- geneSets
    gsc <- gscList
  }
  if (addUserInfo == "yes") {
    tmp <- try(addInfo <- as.data.frame(addInfo, stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("failed to convert additional info in argument 'addInfo' into a data.frame")
    }
  }
  if (class(addInfo) == "data.frame") {
    if (ncol(addInfo) != 2) 
      stop("additional info in argument 'file' or 'addInfo' has to contain 2 columns")
    tmp <- nrow(addInfo)
    addInfo <- unique(addInfo[addInfo[, 1] %in% names(gsc), 
                              ])
  }
  else {
  }
  res <- list(gsc, addInfo)
  names(res) <- c("gsc", "addInfo")
  class(res) <- "GSC"
  return(res)
}



plot_gsa_hyper_heatmap <- function(cds, gsa_results, significance=0.05)
{
  hyper_df <- ldply(gsa_results, function(gsa_res)
  {
    data.frame(gene_set = names(gsa_res$pvalues), pval = gsa_res$pvalues, qval = gsa_res$p.adj)
  })
  colnames(hyper_df)[1] <- "cluster_id"
  
  #hyper_df 
  
  hyper_df <- subset(hyper_df, qval <= significance)
  print (head(hyper_df))
  hyper_df <- merge(hyper_df, ddply(hyper_df, .(gene_set), function(x) { nrow(x) }), by="gene_set")
  #print (hyper_df)
  hyper_df$gene_set <- factor(hyper_df$gene_set, levels=unique(arrange(hyper_df, V1, cluster_id)$gene_set))
  
  qplot(cluster_id, gene_set, fill=-log10(qval), geom="tile", data=hyper_df) + scale_fill_gradientn(colours=rainbow(7))
}



collect_gsa_hyper_results <- function(cds, gsc, clusters)
{
  gene_universe <- unique(as.character(fData(cds)$gene_short_name))
  gsa_results <- list()
  cluster_ids <- unique(clusters)
  for (i in (1:length(cluster_ids))) {
    cluster_genes <- unique(fData(cds[names(clusters[clusters == i]),])$gene_short_name)
    gsaRes <- runGSAhyper(cluster_genes, gsc=gsc, universe=gene_universe)
    gsa_results[[length(gsa_results) + 1]] <- gsaRes
  }
  names(gsa_results) <- cluster_ids
  gsa_results
}

filter_gsc <- function(gsc, cds, alias_table, allowed_genes, gsc_whitelist)
{
  gene_gsc_names_allowed <- sapply(names(gsc$gsc), function(gsc_name, cds, alias_table, allowed_genes) { 
    
    allowed_gene_symbols <- fData(cds[allowed_genes,])$gene_short_name
    
    gsc_matches_to_official_symbol <- sum(grepl(paste("^",gsc_name, sep=""), allowed_gene_symbols))
    gsc_matches_to_alias <- sum(grepl(paste("^", gsc_name, sep=""), subset(alias_table, symbol %in% allowed_gene_symbols)$alias_symbol))
    allowed_gsc_names <- gsc_matches_to_official_symbol + gsc_matches_to_alias > 0
    allowed_gsc_names
  },
  cds, alias_table, allowed_genes
  )
  gsc$gsc <- gsc$gsc[union(names(gene_gsc_names_allowed)[gene_gsc_names_allowed], gsc_whitelist)]
  gsc
}

