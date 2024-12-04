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
# Summary: Identify markers for Visium HD broad cell types in Seurat
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
  library("cowplot")
})

# Define output folder
output_dir <- "/path/to/clustering/output/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load integrated HD object
s <- readRDS("/path/to/integrated/HD/object.rds")

# Remove undefined clusters
s <- subset(s, broad_celltype != "Unknown")

# Recorrect SCT data
DefaultAssay(s) <- "SCT"
s <- PrepSCTFindMarkers(s)

# Set idents
Idents(s) <- s$broad_celltype

# Find positive markers for each cluster
markers <- list()
for (cluster in unique(s$broad_celltype)){
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
write.csv(markers, paste0(output_dir, "markers.csv"), row.names = FALSE, quote = FALSE)

