# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AN1792 Project                              -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Written by: Anne Forsyth
# Summary: Perform hierarchical clustering on lecanemab LOESS trajectories  
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ggrepel")
  library("ggthemes")
  library("ggpubr")
  library("Seurat")
  library("factoextra")
  library("scales")
  library("patchwork")
  library("ggnewscale")
  library("clusterProfiler")
  library("pheatmap")
  library("colorspace")
  library("stringi")
  library("randomcoloR")
  library("UpSetR")
})

# Define output folder
output_folder <- "/path/to/lecanemab/loess/"

# Load and scale predictions
LCMB <- read.csv(paste0(output_folder, "data/LCMB_all_gene_predictions_broad.csv"))
LCMB_scaled <- LCMB
LCMB_scaled[,colnames(LCMB_scaled) != "amyloid"] <- scale(LCMB_scaled[,colnames(LCMB_scaled) != "amyloid"])
LCMB_scaled <- LCMB_scaled[,colSums(is.na(LCMB_scaled)) == 0]

# Perform hierarchical clustering and generate k = 2:12 clusters 
LCMB_dist <- dist(t(LCMB_scaled[,2:ncol(LCMB_scaled)]))
LCMB_clust <- hclust(LCMB_dist, method = "complete")
cluster_nums <- 2:12
LCMB_clust_list <- list()
for (k in cluster_nums) {
  cur_cut <- data.frame(cutree(LCMB_clust, k = k))
  LCMB_clust_list[[k]] <- cur_cut
}

# Save cluster data
LCMB_clust_merged <- as.data.frame(do.call(cbind, LCMB_clust_list[cluster_nums]))
colnames(LCMB_clust_merged) <- paste0("nclust_", cluster_nums)
write.csv(LCMB_clust_merged, paste0(output_folder, "data/LCMB_clusters_all_genes_broad.csv"))

# Load cluster annotations and format data for plotting
LCMB_clust_merged <- read.csv(paste0(output_folder, "data/LCMB_clusters_all_genes_broad.csv"), row.names = 1)
LCMB_long <- pivot_longer(LCMB_scaled, colnames(LCMB_scaled)[2:ncol(LCMB_scaled)], names_to = "gene", values_to = "expression")

# Define uniform y axis limits
min <- min(LCMB_long$expression)
max <- max(LCMB_long$expression)

# Generate plots for select numbers of clusters (11 are ultimately used)
for (n_clust in c(6, 10, 11, 12)) {
  
  # Define nrow and ncol for ggarrange
  if (n_clust == 6) {
    nrow <- 2
    ncol <- 3
  } else if (n_clust == 10) {
    nrow <- 2
    ncol <- 5
  } else {
    nrow <- 3
    ncol <- 4
  }
  
  # Define cluster colors
  colors <- hue_pal()(n_clust)
  
  # Generate line plots for each cluster, with an additional LOESS curve for average predicted expression
  plt_list <- list()
  for (clust in 1:n_clust) {
    cur_data <- LCMB_long[LCMB_long$gene %in% row.names(LCMB_clust_merged)[LCMB_clust_merged[,paste0("nclust_", n_clust)] == clust],]
    
    if (length(unique(cur_data$gene)) > 1) {
      avg_data <- cur_data %>% group_by(amyloid) %>% summarize(average = mean(expression))
      
      plt <- ggplot(cur_data, aes(x = amyloid, y = expression)) +
        theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black')) +
        geom_line(linewidth = 0.3, aes(group = gene), color = colors[clust]) +
        geom_line(data = avg_data,
                  aes(x = amyloid, y = average),
                  stat="smooth", method = "loess", span = 0.75, se = FALSE,
                  color = darken(colors[[clust]], amount = 0.3), linewidth = 2) +
        theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
              axis.text.x = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
              axis.text.y = element_text(size = 15)) +
        ggtitle(paste0("Cluster ", clust, ": ", length(unique(cur_data$gene)), " Genes")) + theme(aspect.ratio = 1) + ylim(min, max) +
        geom_vline(xintercept = 183, linetype = 3)
    } else {
      plt <- ggplot(cur_data, aes(x = amyloid, y = expression)) +
        theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black')) +
        geom_line(linewidth = 0.3, aes(group = gene), color = colors[clust]) +
        theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
              axis.text.x = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
              axis.text.y = element_text(size = 15)) +
        ggtitle(paste0("Cluster ", clust, ": ", length(unique(cur_data$gene)), " Genes")) + theme(aspect.ratio = 1) + ylim(min, max) +
        geom_vline(xintercept = 183, linetype = 3)
    }
    plt_list[[clust]] <- plt
  }
  
  # Save plot grid
  pdf(paste0(output_folder, "plots/LCMB_clusters_", n_clust, "_all_genes.pdf"), height = 10*nrow, width =10*ncol)
  print(ggarrange(plotlist = plt_list, ncol = ncol, nrow = nrow))
  dev.off()
}











