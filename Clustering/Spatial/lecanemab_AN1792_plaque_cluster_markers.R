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
# Summary: Identify markers for lecanemab and AN1792 Aß clusters
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
})

# Define clustering resolution 
res <- 0.4

# Define output folder
output_dir <- paste0("/path/to/clustering/output/harmony_", res, "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load all-cohorts gray matter amyloid-rich integrated Seurat object and recorrect SCT data
s <- readRDS("/path/to/integrated/all/cohorts/amyloid/object.rds")
DefaultAssay(s) <- "SCT"
for (name in names(s@assays$SCT@SCTModel.list)) {
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}
s <- PrepSCTFindMarkers(s)

# Create temporary variable for clusters and set idents
s$clusters <- s@meta.data[[paste0("harmony_snn_res.", res)]]
Idents(s) <- "clusters"

# Find positive markers for each cluster 
markers <- list()
for (cluster in unique(s$clusters)){
  
  # Test genes expressed in 25% of cluster
  LFC <- FoldChange(s, ident.1 = cluster, slot = "data", assay = "SCT", base = 2)
  genes <- rownames(LFC)[LFC$pct.1 > 0.25]
  
  # Exclude contamination genes
  if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes)) != 0) {
    genes <- genes[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes)] 
  }
  
  # Find positive cluster markers with significant adjusted p value
  results <- FindMarkers(s, ident.1 = cluster, features = genes, recorrect_umi = FALSE, min.pct = 0, only.pos = TRUE, 
                         fc.slot = "data", base = 2)
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
}

# Compile and save results
all_markers <- data.table::rbindlist(markers)|>as.data.frame()
write.csv(all_markers, paste0(output_dir, "markers.csv"), row.names = FALSE, quote = FALSE)

