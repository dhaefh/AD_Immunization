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
# Summary: Doublet detection for lecanemab with DoubletFinder
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("ggthemes")
  library("ggrepel")
  library("grid")
  library("DoubletFinder")
  library("doMC")
  library("xlsx")
  library("RColorBrewer")
})

# Define output folder
output_dir <- "/path/to/lecanemab/QC/output/"

# Define plot folder
plot_dir <- paste0(output_dir, "plots/DoubletFinder/")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Define sample object subfolder 
object_dir <- paste0(output_dir, "data/DoubletFinder/")
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

# Load Seurat object with low quality cells removed
s <- readRDS(paste0(output_dir, "data/s_lowquality_removed.rds"))

# Define estimated multiplet rates per pool
doublet_formation_rates <- c(0.064, 0.064, 0.052)
names(doublet_formation_rates) <- c("pool1", "pool2", "CAA")

# Define function to run DoubletFinder
run_doubletfinder <- function(s) {
  
  # Get pool-specific doublet formation rate
  pool <- unique(s$pool)
  doublet_formation_rate <- doublet_formation_rates[pool]
  print(paste0(pool, " rate: ", doublet_formation_rate))
  
  # Standard normalization and scaling 
  s <- s %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
  
  # Default Seurat clustering and TSNE 
  s <- s %>% RunPCA() %>% FindNeighbors(dims = 1:10) %>% FindClusters()
  s <- RunTSNE(s, dims = 1:10, check_duplicates = FALSE) 
  
  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(s, PCs = 1:10, sct = FALSE) 
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  max_index <- which.max(bcmvn$BCmetric)
  optimal_pK <- as.numeric(as.character(bcmvn[max_index, "pK"]))
  
  # Homotypic Doublet Proportion Estimate
  annotations <- s@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_formation_rate*nrow(s@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # Run DoubletFinder without homotypic adjustment 
  s <- doubletFinder(s, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # Rename meta column
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data))] <- "DF.unadj"
  
  # Run DoubletFinder with homotypic adjustment
  pANN <- colnames(s@meta.data)[grep("^pANN", colnames(s@meta.data))]
  print(pANN)
  s <- doubletFinder(s, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi.adj, reuse.pANN = pANN, sct = FALSE)
  
  # Rename meta column
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data))] <- "DF.adj"
  
  # Plot unadjusted vs. adjusted doublets in TSNE coordinates
  plt <- (DimPlot(s, reduction = "tsne", group.by = "DF.unadj", pt.size = 1, shuffle = TRUE, alpha = 0.6) + 
            scale_color_manual(values = list(Doublet = "darkblue", Singlet = "lightgray")) + theme(aspect.ratio = 1) + ggtitle("Unadjusted")) + 
    (DimPlot(s, reduction = "tsne", group.by = "DF.adj", pt.size = 1, shuffle = TRUE, alpha = 0.6) + 
       scale_color_manual(values = list(Doublet = "darkblue", Singlet = "lightgray")) + theme(aspect.ratio = 1) + ggtitle("Adjusted for Homotypic Proportion"))
  
  pdf(paste0(plot_dir, unique(s$sample_id), "_tsne.pdf"), height = 10, width = 20)
  print(plt)
  dev.off()
  
  return(s)
}

# Create list of sample objects
s_list <- SplitObject(s, split.by = "sample_id")

# Run DoubletFinder on each sample
for (s in s_list) {
  s <- run_doubletfinder(s)
  saveRDS(s, paste0(object_dir, unique(s$sample_id), ".rds"))
}

