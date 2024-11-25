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
# Summary: Identify positive and negative markers for Aß cluster 6
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("DESeq2")
  library("ggrepel")
  library("UpSetR")
  library("randomcoloR")
  library("ComplexUpset")
})

# Define output folder
output_dir <- "/path/to/cluster6/marker/output/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load all-cohorts gray matter Aß-rich integrated Seurat object and recorrect SCT data
s <- readRDS("/path/to/integrated/all/cohorts/amyloid/object.rds")
DefaultAssay(s) <- "SCT"
for (name in names(s@assays$SCT@SCTModel.list)) {
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}
s <- PrepSCTFindMarkers(s)

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Set idents
Idents(s) <- "harmony_snn_res.0.4"

# Exclude contamination genes from testing 
genes <- rownames(s)
if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes)) != 0) {
  genes <- genes[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes)] 
}

# Identify positive and negative markers for Aß cluster 6 vs. all other Aß clusters
results <- FindMarkers(s, ident.1 = 6, min.pct = 0.1, recorrect_umi = FALSE, only.pos = FALSE, features = genes, 
                       fc.slot = "data", base = 2, logfc.threshold = -Inf)
results$BH <- p.adjust(results$p_val, method = "BH")
results$gene <- rownames(results)

# Define DE markers 
results$DE <- "Not DE"
results$DE[results$avg_log2FC > fc_thresh & results$BH < p_thresh] <- "Upregulated"
results$DE[results$avg_log2FC < -fc_thresh & results$BH < p_thresh] <- "Downregulated"
results$DE_gene <- NA
results$DE_gene[results$DE != "Not DE"] <- results$gene[results$DE != "Not DE"]

# Save results
write.csv(results, paste0(output_dir, "markers.csv"))


