require(monocle3)
library(tidyr)
library(tidymodels)

test_drug_dependence <- function(cds_subset, cds, reference_cells=NULL, min_fraction_cells=0.01, pseudodose=0.01, residualModelFormulaTerms=NULL){
  #print(cell_ids_to_test)
  cell_ids_to_test = as.character(cds_subset$cell)
  
  print(length(cell_ids_to_test))
  
  cds_subset = cds[,reference_cells]
  cds_subset = detect_genes(cds_subset)
  genes_in_reference_set = row.names(subset(fData(cds_subset), num_cells_expressed > ncol(cds_subset)*min_fraction_cells))
  
  #print(length(genes_in_reference_set))
  #print (head(genes_in_reference_set))
  cds_subset = cds[,cell_ids_to_test]
  cds_subset = detect_genes(cds_subset)
  genes_in_treated_set = row.names(subset(fData(cds_subset), num_cells_expressed > ncol(cds_subset)*min_fraction_cells))
  print(length(genes_in_treated_set))
  #print (head(genes_in_treated_set))
  cds_subset = cds[,union(cell_ids_to_test, reference_cells)]
  #print(cds_subset)
  #cds_subset = estimateSizeFactors(cds_subset)
  #cds_subset = estimateDispersions(cds_subset)
  cds_subset = cds_subset[union(genes_in_reference_set, genes_in_treated_set),]
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
  
  #detect_genes(cds)
  sci_chem_model_tbl = fit_models(cds_subset, model_formula_str = modelFormulaStr, reltol=1e-5, verbose=TRUE)
  sci_chem_model_coefs = coefficient_table(sci_chem_model_tbl)
  sci_chem_model_coefs
  #data.frame(xxx=1)
}
