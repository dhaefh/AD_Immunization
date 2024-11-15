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
# Summary: Cluster Visium HD data and find cluster markers in Seurat
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
})

# Define output folder
output_folder <- "/path/to/clustering/output/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Load integrated HD object
s <- readRDS("/path/to/integrated/HD/object.rds")

# Find clusters for a sequence of resolutions
s <- FindNeighbors(s, reduction = "harmony", dims = 1:30, graph.name = c('harmony_nn', 'harmony_snn'))
for (i in seq(0.05, 0.5, 0.05)) {
  s <- FindClusters(s, resolution = i, graph.name = 'harmony_snn')
}

# Save object
saveRDS(s, "/path/to/integrated/HD/object.rds")

# Plot clusters for each resolution
plt_list <- list()
index <- 1
for (i in seq(0.05, 0.5, 0.05)) {
  colname <- paste0("harmony_snn_res.", i)
  s$temp <- factor(s@meta.data[[colname]], levels = 0:(length(unique(s@meta.data[[colname]])) - 1))
  
  plt <- DimPlot(s, group.by = "temp", reduction = "harmonyumap", raster = FALSE, shuffle = TRUE) + 
    scale_color_manual(values = randomcoloR::distinctColorPalette(length(unique(s$temp)))) +
    ggtitle(paste0("Resolution = ", i)) + theme(aspect.ratio = 1)
  
  plt_list[[index]] <- plt
  
  index <- index + 1
}
pdf(paste0(output_folder, "harmony_cluster_umaps.pdf"), height = 20, width = 50)
ggpubr::ggarrange(plotlist = plt_list, nrow = 2, ncol = 5)
dev.off()

# Create subfolder for FindMarkers output 
marker_dir <- paste0(output_folder, "harmony_markers_0.2/")
dir.create(marker_dir, showWarnings = FALSE, recursive = TRUE)

# Recorrect SCT data
DefaultAssay(s) <- "SCT"
s <- PrepSCTFindMarkers(s)

# Create temporary variable for clusters and set idents
s$clusters <- s$harmony_snn_res.0.2
Idents(s) <- s$clusters

# Find positive markers for each cluster  
markers <- list()
for (cluster in unique(s$clusters)){
  print(cluster)
  
  tryCatch(
    {
      # Test genes expressed in 1% of cluster
      LFC <- FoldChange(s, ident.1 = cluster, slot = "data", assay = "SCT", base = 2)
      genes <- rownames(LFC)[LFC$pct.1 > 0.01]
      
      # Exclude contamination genes
      if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes)) != 0) {
        genes <- genes[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes)]
      }
      
      # Find positive cluster markers with significant adjusted p value
      results <- FindMarkers(s, ident.1 = cluster, features = genes, recorrect_umi = FALSE, min.pct = 0, only.pos = TRUE)
      results$BH <- p.adjust(results$p_val, method = "BH")
      results$cluster <- cluster
      results$gene <- rownames(results)
      results <- dplyr::filter(results, BH < 0.05)
      
      # Ensure PFC is finite; if all adjusted p values are zero set neg_logBH to a constant so genes can be ranked by PFC
      results$neg_logBH <- -log10(results$BH)
      if (sum(results$neg_logBH != Inf) == 0) {
        results$neg_logBH <- 500
      } else {
        results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
      }
      results$PFC <- results$neg_logBH*abs(results$avg_log2FC)
      markers[[cluster]] <- results
    },
    error = function(e) {
      print(paste0("Error in cluster ", cluster))
    }
  )
}

# Compile and save results
markers <- data.table::rbindlist(markers)|>as.data.frame()
write.csv(markers, paste0(marker_dir, "markers.csv"), row.names = FALSE, quote = FALSE)

