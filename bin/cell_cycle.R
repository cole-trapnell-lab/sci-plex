estimate_cell_cycle = function(cds, g1s_markers, g2m_markers){
  cds_g1s = cds[fData(cds)$gene_short_name %in% g1s_markers,] 
  aggregate_g1s_expression = exprs(cds_g1s)
  aggregate_g1s_expression = t(t(aggregate_g1s_expression) / pData(cds_g1s)$Size_Factor)
  aggregate_g1s_expression = Matrix::colSums(aggregate_g1s_expression)
  
  cds_g2m = cds[fData(cds)$gene_short_name %in% g2m_markers,] 
  aggregate_g2m_expression = exprs(cds_g2m)
  aggregate_g2m_expression = t(t(aggregate_g2m_expression) / pData(cds_g2m)$Size_Factor)
  aggregate_g2m_expression = Matrix::colSums(aggregate_g2m_expression)
  #pData(cds)$cell_cycle_score = aggregate_g1s_expression
  pData(cds)$g1s_score = log(aggregate_g1s_expression+1)
  pData(cds)$g2m_score = log(aggregate_g2m_expression+1)
  pData(cds)$proliferation_index = log(aggregate_g1s_expression + aggregate_g2m_expression + 1)
  return(cds)
}


aggregate_gene_score = function(cds, cycle_markers=c("CCNB1", "CCNB2","CDK1")){
  cds_cell_cycle_norm = cds[fData(cds)$gene_short_name %in% cycle_markers,] 
  aggregate_marker_expression = exprs(cds_cell_cycle_norm)
  aggregate_marker_expression = t(t(aggregate_marker_expression) / pData(cds_cell_cycle_norm)$Size_Factor)
  aggregate_marker_expression = log10(aggregate_marker_expression + 1)
  aggregate_marker_expression = t(scale(t(aggregate_marker_expression)))
  aggregate_marker_expression = aggregate_marker_expression[rowSums(is.na(aggregate_marker_expression)) == 0,]
  aggregate_marker_expression = colMeans(aggregate_marker_expression)
  pData(cds)$cell_cycle_score = aggregate_marker_expression
  return(cds)
}