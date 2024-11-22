# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                         AN1792 Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Written by: Anne Forsyth 
# Summary: Generate heatmap of marker gene expression
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
})

# Load Seurat object 
s <- readRDS("/path/to/Seurat/object.rds")

# Recorrect SCT data
if ("Spatial" %in% names(s@assays)) {
  for (name in names(s@assays$SCT@SCTModel.list)) {
    s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
  }
}
DefaultAssay(s) <- "SCT"
s <- PrepSCTFindMarkers(s)

# Define cluster variable
s$cluster <- s@meta.data[["<name of cluster column>"]]

# Load FindMarkers output 
markers <- read.csv("/path/to/markers.csv")

# Define marker genes to plot
marker_genes_to_plot <- c("<list of marker genes>")

# Optionally, select top n markers per cluster based on probability fold change 
n <- 2
marker_genes_to_plot <- c()
for(clust in all_markers$cluster|>unique()){
  data <- dplyr::filter(all_markers, cluster==clust)
  print(sum(data$PFC == Inf))
  data <- data %>% dplyr::arrange(desc(PFC))
  marker_genes_to_plot <- c(marker_genes_to_plot, data$gene|>head(n))
}
marker_genes_to_plot <- marker_genes_to_plot %>% unique()

# Calculate average recorrected SCT expression per cluster
data <- GetAssayData(s, assay = "SCT", layer = "data")
data <- data[marker_genes_to_plot,]
data <- as.matrix(data)
colnames(data) <- s$cluster
aggregated <- t(rowsum(t(data), group = colnames(data)))
for (col in colnames(aggregated)) {
  aggregated[,col] <- aggregated[,col] / sum(colnames(data) == col)
}

# Standardize expression across clusters
scaled <- t(aggregated)
scaled <- scale(scaled)
scaled <- t(scaled)

# Format data for plotting
data <- as.data.frame(scaled) %>% rownames_to_column(var = "gene") %>% tidyr::gather(key = "cluster",value = "expression",-gene)

# Cluster rows (genes)
row_clust <- hclust(dist(as.matrix(scaled)))
roworder <- row_clust$labels[row_clust$order]
data$gene <- factor(data$gene, levels = roworder)

# Cluster columns (clusters)
col_clust <- hclust(dist(as.matrix(t(scaled))))
colorder <- col_clust$labels[col_clust$order]
data$cluster <- factor(data$cluster, levels = colorder)

# Generate continuous x and y variables  
data$y <- NA
i <- 1
for (gene in levels(data$gene)) {
  data$y[data$gene == gene] <- i
  i <- i + 1
}

data$x <- NA
i <- 1
for (cluster in levels(data$cluster)) {
  data$x[data$cluster == cluster] <- i
  i <- i + 1
}

# Initialize heatmap
plt <- ggplot(data, mapping = aes(x = x, y = y)) + geom_tile(aes(fill = expression), color = "white", linewidth = 0.5) +
  scale_fill_gradientn(colours=c("#330099", "#99CCFF", "#FFFFCC", "#FFCC66", "#CC0000"))

# Add dendrogram corresponding to genes - adjust xlim arguments based on number of clusters on x axis
plt <- plt + ggdendroplot::geom_dendro(row_clust, pointing = "side", xlim=c(9.5, 10.5))

# Add dendrogram corresponding to clusters - adjust ylim arguments based on number of genes on y axis
plt <- plt + ggdendroplot::geom_dendro(col_clust, pointing = "updown", ylim=c(31.5, 32.25))  

# Add formatting
plt <- plt + theme(panel.background = element_rect(fill = 'white')) +
  ggtitle("<plot title>") + # Specify title
  theme(plot.title = element_text(hjust = 0.5, vjust = 2), axis.title = element_blank(), 
        axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust=1)) + 
  coord_equal() + labs(fill = "Scaled Expression")



