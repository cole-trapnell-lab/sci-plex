## This function was swiped from DESeq (Anders and Huber) and modified for our purposes
#' @importFrom stats glm Gamma
parametricDispersionFit <- function( disp_table, verbose = FALSE, initial_coefs=c(1e-6, 1) )
{
  coefs <- initial_coefs
  iter <- 0
  while(TRUE) {
    residuals <- disp_table$disp / ( coefs[1] + coefs[2] / disp_table$mu )
    good <- disp_table[which(disp_table$disp > 0 & (residuals < 10000) ),]
    #good <- disp_table
    if(verbose)
      fit <- glm( disp ~ I(1/mu), data=good,
                  family=Gamma(link="identity"), start=coefs )
    else
      suppressWarnings(fit <- glm( disp ~ I(1/mu), data=good,
                                   family=Gamma(link="identity"), start=coefs ))
    
    oldcoefs <- coefs
    coefs <- coefficients(fit)
    if (coefs[1] < initial_coefs[1]){
      coefs[1] <- initial_coefs[1]
    }
    if (coefs[2] < 0){
      stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')" )
      # stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation with DESeq package. \nYou may try to use the following code for the calculation: 
      #     cds_deseq <- DESeq::newCountDataSet(round(exprs(cds)), pData(cds)$Hours, sizeFactors = sizeFactors(cds))\n 
      #     cds_deseq <- estimateDispersions(cds, fitType = 'local')\n
      #     cds@dispFitInfo$blind <- monocle_HSMM_deseq@fitInfo$pooled" )
    }
    #     if( !all( coefs > 0 ) ){
    #       #print(data.frame(means,disps))
    #       stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')" )
    #     }
    if( sum( log( coefs / oldcoefs )^2 ) < coefs[1] )
      break
    iter <- iter + 1
    #print(coefs)
    if( iter > 10 ) {
      warning( "Dispersion fit did not converge." )
      break
    }
  }
  
  if( !all( coefs > 0 ) ){
    #print(data.frame(means,disps))
    stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')" )
  }
  
  #names( coefs ) <- c( "asymptDisp", "extraPois" )
  #ans <- function( q )
  #  coefs[1] + coefs[2] / q
  #ans
  #coefs
  #rm(ans)
  rm(disp_table)
  #rm(cds_pdata)
  #rm(res)
  #rm(cds)
  
  list(fit, coefs)
}


#' @importFrom Biobase exprs pData fData
#' @importFrom DelayedArray DelayedArray
disp_calc_helper_NB <- function(cds, expressionFamily, min_cells_detected){
  
  
  rounded <- round(exprs(cds))
  nzGenes <- Matrix::rowSums(rounded > 0.1) #cds@lowerDetectionLimit)
  nzGenes <- names(nzGenes[nzGenes > min_cells_detected])
  
  # Note: we do these operations as DelayedArray ops because standard operations will trigger a conversion 
  # to an in-memory dense matrix. DelayedArray uses block processing. Control block size with:
  # options(DelayedArray.block.size=100e6) 
  # We should make this clear in the documentation, and possibly
  # emit a message to users on calling this (and possibly other) functions.
  
  #Progress Bar here
  x <- DelayedArray(t(t(rounded[nzGenes,]) / size_factors(cds[nzGenes,])))
  
  xim <- mean(1/ size_factors(cds[nzGenes,]))
  

    f_expression_mean <- as(DelayedMatrixStats::rowMeans2(x), "sparseVector")

  
  # For NB: Var(Y)=mu*(1+mu/k)
  f_expression_var <- DelayedMatrixStats::rowVars(x)
  
  disp_guess_meth_moments <- f_expression_var - xim * f_expression_mean
  
  disp_guess_meth_moments <- disp_guess_meth_moments / (f_expression_mean^2) #fix the calculation of k
  
  
  res <- data.frame(mu=as.vector(f_expression_mean), disp=as.vector(disp_guess_meth_moments))
  res[res$mu == 0]$mu = NA
  res[res$mu == 0]$disp = NA
  res$disp[res$disp < 0] <- 0
  
  res <- cbind(gene_id=row.names(fData(cds[nzGenes,])), res)
  res
}


#' Helper function to estimate dispersions
#' @importFrom Biobase pData
#' @importFrom stats cooks.distance
#' @importFrom stringr str_split str_trim
#' @importFrom dplyr %>%
#' @param cds a CellDataSet that contains all cells user wants evaluated
#' @param min_cells_detected Only include genes detected above lowerDetectionLimit in at least this many cells in the dispersion calculation
#' @param removeOutliers a boolean it determines whether or not outliers from the data should be removed
#' @param verbose Whether to show detailed running information.
estimateDispersionsForCellDataSet <- function(cds, min_cells_detected, removeOutliers, verbose = FALSE)
{
  mu <- NA
  
  progress_opts <- options()$dplyr.show_progress
  options(dplyr.show_progress = T)
  
  cds_pdata <- dplyr::group_by_(dplyr::select_(rownames_to_column(as.data.frame(pData(cds))), "rowname"))
  disp_table <- as.data.frame(disp_calc_helper_NB(cds, cds@expressionFamily, min_cells_detected))
  #disp_table <- data.frame(rowname = row.names(type_res), CellType = type_res)
  
  #message("fitting disersion curves")
  #print (disp_table)
  if(!is.list(disp_table))
    stop("Parametric dispersion fitting failed, please set a different lowerDetectionLimit")
  #disp_table <- do.call(rbind.data.frame, disp_table)
  disp_table <- subset(disp_table, is.na(mu) == FALSE)
  res <- parametricDispersionFit(disp_table, verbose)
  fit <- res[[1]]
  coefs <- res[[2]]
  #removeOutliers = TRUE
  if (removeOutliers){
    CD <- cooks.distance(fit)
    #cooksCutoff <- qf(.99, 2, ncol(cds) - 2)
    cooksCutoff <- 4/nrow(disp_table)
    #print (head(CD[CD > cooksCutoff]))
    #print (head(names(CD[CD > cooksCutoff])))
    message (paste("Removing", length(CD[CD > cooksCutoff]), "outliers"))
    outliers <- union(names(CD[CD > cooksCutoff]), setdiff(row.names(disp_table), names(CD)))
    res <- parametricDispersionFit(disp_table[row.names(disp_table) %in% outliers == FALSE,], verbose)
    fit <- res[[1]]
    coefs <- res[[2]]
  }
  names( coefs ) <- c( "asymptDisp", "extraPois" )
  ans <- function( q )
    coefs[1] + coefs[2] / q
  attr( ans, "coefficients" ) <- coefs
  
  res <- list(disp_table = disp_table, disp_func = ans)
  return(res)
}

dispersionTable <- function (disp_results) 
{
  if (is.null(disp_results)) {
    warning("Warning: estimateDispersions only works, and is only needed, when you're using a CellDataSet with a negbinomial or negbinomial.size expression family")
    stop("Error: no dispersion model found. Please call estimateDispersions() before calling this function")
  }
  disp_df <- data.frame(gene_id = disp_results$disp_table$gene_id, 
                        mean_expression = disp_results$disp_table$mu, 
                        dispersion_fit = disp_results$disp_func(disp_results$disp_table$mu), 
                        dispersion_empirical = disp_results$disp_table$disp)
  return(disp_df)
}
