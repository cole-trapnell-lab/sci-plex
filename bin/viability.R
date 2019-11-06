# Estimate viability:

estimate_viability_helper = function(pData_rows, sci_chem_models){
  # Make sure we are estimating viability for a single cell type
  ct = unique(pData_rows$cell_type)
  stopifnot(length(ct) == 1)
  ct = ct[1]
  
  # Make sure we are estimating viability for a single treatment
  tmt = unique(pData_rows$treatment)
  stopifnot(length(tmt) == 1)
  tmt = tmt[1]
  
    
  sci_chem_model = 
    sci_chem_models %>% 
    filter(cell_type == ct & treatment == tmt) %>% 
    dplyr::select(drc_model)
  
  
  stopifnot(length(sci_chem_model) == 1)
  
  res = tryCatch({
    #print(sci_chem_model)
    sci_chem_model = sci_chem_model$drc_model[[1]]
    #print (head(pData_rows))
    newdata =  data.frame(pData_rows$dose)
    
    #newdata = data.frame(dose=seq(min(model$origData$dose), max(model$origData$dose),length.out=100))
    cell_counts = predict(sci_chem_model, newdata=newdata)
    #eps_dose = min(model$origData$dose[model$origData$dose>0], na.rm=T)
    eps_dose = min(sci_chem_model$origData$dose[sci_chem_model$origData$dose>0], na.rm=T)
    newdata =  data.frame(dose=exp(log(min(sci_chem_model$origData$dose, na.rm=T) + eps_dose/2)))
    max_cell_counts = predict(sci_chem_model, newdata=newdata)[1]
    viability = cell_counts / max_cell_counts
    names(viability) = pData_rows$cell
    viability
  }, error = function(e) {
    print(e)
    viability = rep(NA, nrow(pData_rows))
    names(viability) = pData_rows$cell
    viability
  })
  
  data.frame(cell = names(res), viability = res)
}

estimate_viability = function(cds, sci_chem_models) {
  viability_df = 
    colData(cds) %>% 
    as.data.frame() %>%
    group_by(cell_type, sample, treatment) %>%
    do(estimate_viability_helper(., sci_chem_models)) 
  
  row.names(viability_df) = viability_df$cell
  
  #print(head(viability_df)) 
  #pData(cds)$viability = viability_df[pData(cds)$cell,]$viability
  
  pData(cds)$viability = pmin(1, viability_df[pData(cds)$cell,]$viability)
  return(cds)
}