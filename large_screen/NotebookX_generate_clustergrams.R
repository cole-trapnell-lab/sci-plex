##############################
# Generate circular clustergrams of the drug relationships from the heatmaps above
##############################

# Create an annotation data.frame
color_df <- drug_annotations %>% 
filter(product_name != "Vehicle" & !is.na(pathway_level_1)) %>% 
dplyr::select(product_name,pathway_level_1)

color_df$product_name <- as.character(color_df$product_name)
color_df$pathway_level_1 <- as.character(color_df$pathway_level_1)


color_df$col <- sapply(color_df$pathway_level_1, function(x){

  if(x == "Antioxidant")return("aquamarine3")
  if(x == "Apoptotic regulation")return("lemonchiffon4")
  if(x == "Cell cycle regulation")return("deepskyblue2")
  if(x == "DNA damage & DNA repair")return("slategray4")
  if(x == "Epigenetic regulation")return("navy")
  if(x == "Focal adhesion signaling")return("brown4")
  if(x == "HIF signaling")return("darkgreen")
  if(x == "JAK/STAT signaling")return("orangered2")
  if(x == "Metabolic regulation")return("gold2")
  if(x == "Neuronal signaling")return("darkolivegreen")
  if(x == "Nuclear receptor signaling")return("chartreuse3")
  if(x == "PKC signaling")return("plum4")
  if(x == "Protein folding & Protein degradation")return("darkorchid2")
  if(x == "TGF/BMP signaling")return("darkcyan")
  if(x == "Tyrosine kinase signaling")return("firebrick1")
  if(x == "Other")return("darkgoldenrod4")
  return(NA)
  })

color_df$col <- as.character(color_df$col)

ann_colors = list("Pathway" = c("Antioxidant"="aquamarine3",
                                "Apoptotic regulation"="lemonchiffon4",
                                "Cell cycle regulation"="deepskyblue2",
                                "DNA damage & DNA repair"="slategray4",
                                "Epigenetic regulation"="navy",
                                "Focal adhesion signaling"="brown4",
                                "HIF signaling"="darkgreen",
                                "JAK/STAT signaling"="orangered2",
                                "Metabolic regulation" = "gold2",
                                "Neuronal signaling" = "darkolivegreen",
                                "Nuclear receptor signaling" = "chartreuse3",
                                "PKC signaling" = "plum4",
                                "Protein folding & Protein degradation" = "darkorchid2",
                                "TGF/BMP signaling" = "darkcyan",
                                "Tyrosine kinase signaling" = "firebrick1",
                                "Other" ="darkgoldenrod4"))

# Isolate dendrogram for A549 and make a data.frame containing annotations
A549_bounded_den <- dendextend::as.ggdend(as.dendrogram(A549_bounded_ph$tree_row))

A549_bounded_den$labels$col <- sapply(A549_bounded_den$labels$label, function(x){
  
  y = as.character(x)
  color <- color_df[color_df$product_name == y,]$col
  return(color)

})

new_label <- sapply(A549_bounded_den$labels$label, function(x){
  
  y = as.character(x)
  new_label <- strsplit(y, "\\s+")[[1]][1]
  return(new_label)

})

A549_bounded_den$labels$new_label <- new_label

# Plot clustergram
ggplot(A549_bounded_den, labels = FALSE) +
scale_y_reverse(expand = c(0.5, 0.2)) +
coord_polar(theta="x") +
geom_text(data = A549_bounded_den$labels, 
              aes(x = x, y = y, label = new_label, 
                color = col, 
                angle = -90 - 360 / length(unique(A549_bounded_den$labels$label)) * seq_along(A549_bounded_den$labels$label)), 
              size = 5.5, hjust = 1) +
ggsave("A549_drug_correlation_matrix_clustergram.png", height = 15, width = 15, dpi = 600)

# Isolate dendrogram for K562 and make a data.frame containing annotations
K562_bounded_den <- dendextend::as.ggdend(as.dendrogram(K562_bounded_ph$tree_row))

K562_bounded_den$labels$col <- sapply(K562_bounded_den$labels$label, function(x){
  
  y = as.character(x)
  color <- color_df[color_df$product_name == y,]$col
  return(color)

})

new_label <- sapply(K562_bounded_den$labels$label, function(x){
  
  y = as.character(x)
  new_label <- strsplit(y, "\\s+")[[1]][1]
  return(new_label)

})

K562_bounded_den$labels$new_label <- new_label

# Plot clustergram
ggplot(K562_bounded_den, labels = FALSE) +
scale_y_reverse(expand = c(0.5, 0.2)) +
coord_polar(theta="x") +
geom_text(data = K562_bounded_den$labels, 
              aes(x = x, y = y, label = new_label, 
                color = col, 
                angle = -90 - 360 / length(unique(K562_bounded_den$labels$label)) * seq_along(K562_bounded_den$labels$label)), 
              size = 5.5, hjust = 1) +
ggsave("K562_drug_correlation_matrix_clustergram.png", height = 15, width = 15, dpi = 600)

# Isolate dendrogram for MCF7 and make a data.frame containing annotations
MCF7_bounded_den <- dendextend::as.ggdend(as.dendrogram(MCF7_bounded_ph$tree_row))

MCF7_bounded_den$labels$col <- sapply(MCF7_bounded_den$labels$label, function(x){
  
  y = as.character(x)
  color <- color_df[color_df$product_name == y,]$col
  return(color)

})

new_label <- sapply(MCF7_bounded_den$labels$label, function(x){
  
  y = as.character(x)
  new_label <- strsplit(y, "\\s+")[[1]][1]
  return(new_label)

})

MCF7_bounded_den$labels$new_label <- new_label

# Plot clustergram
ggplot(MCF7_bounded_den, labels = FALSE) +
scale_y_reverse(expand = c(0.5, 0.2)) +
coord_polar(theta="x") +
geom_text(data = MCF7_bounded_den$labels, 
              aes(x = x, y = y, label = new_label, 
                color = col, 
                angle = -90 - 360 / length(unique(MCF7_bounded_den$labels$label)) * seq_along(MCF7_bounded_den$labels$label)), 
              size = 5.5, hjust = 1) +
ggsave("MCF7_drug_correlation_matrix_clustergram.png", height = 15, width = 15, dpi = 600)