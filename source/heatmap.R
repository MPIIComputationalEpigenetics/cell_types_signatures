library(gplots)
library(RColorBrewer)

plot_heatmap <- function(signatures, mean_df_matrix, experiments_meta, 
                         signature_score_numeric, filename)
{
  signature_regions <- rbindlist(signatures) %>% select(CHROMOSOME, START, END) %>% distinct()
  
  
  # filter the DNA methylation matrix for those signature regions
  signature_score_matrix <- dplyr::semi_join(mean_df_matrix, 
                                             signature_regions, 
                                             by = c("CHROMOSOME", "START", "END")) 
  signature_score_numeric <- as.matrix(signature_score_matrix[,4:ncol(signature_score_matrix)])
  signature_score_numeric[is.na(signature_score_numeric)] <- 0
  
  rownames(experiments_meta) <- experiments_meta$name
  unique_biosources <- unique(experiments_meta$user_celltype)

  # define cell types and colors
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  
  color_map <- data.frame(biosource = unique_biosources,
                          color = getPalette(length(unique_biosources)))
  
  exp_names <- colnames(signature_score_numeric)
  biosource_colors <- data.frame(name = exp_names, biosource = experiments_meta[exp_names, "user_celltype"])
  biosource_colors <- left_join(biosource_colors, color_map, by = "biosource")
  color_vector <- as.character(biosource_colors$color)
  names(color_vector) <-  biosource_colors$biosource
  
  #plot heatmap using pearson correlation for hierarchical clustering
  pdf(file = paste0(filename, "_heatmap.pdf"), paper = "a4")
  heatmap.2(signature_score_numeric,labRow = NA, labCol = NA,
            trace = "none", ColSideColors = color_vector,
            hclust=function(x) hclust(x,method="complete"),
            distfun=function(x) as.dist(1-cor(t(x), method = "pearson")), Rowv = TRUE, dendrogram = "column",
            key.xlab = "beta value", denscol = "black", keysize = 1.5,
            key.par = list(mar = c(8.5, 2.5, 1, 1)), key.title = NA)
  
  # Next, we add a legend showing which cell type has which color
  plot.new()
  legend(x = 0, y = 1,
         legend = color_map$biosource,
         col = as.character(color_map$color),
         text.width = 0.6,
         lty= 1,
         lwd = 6,
         cex = 0.7,
         y.intersp = 0.7,
         x.intersp = 0.7,
         inset=c(-0.21,-0.11))
  dev.off()
}
