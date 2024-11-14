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
# Summary: Initial clustering of scRNAseq data in Seurat
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
output_dir <- "/path/to/clustering/output/"

# Create plot output folder
plot_dir <- paste0(output_dir, "plots/harmony/")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load integrated scRNAseq Seurat object
s <- readRDS("/path/to/integrated/scRNAseq/object.rds")

# Find clusters for a sequence of resolutions
s <- FindNeighbors(s, reduction = "harmony", dims = 1:30, graph.name = c('harmony_nn', 'harmony_snn'))
for (i in seq(0.25, 1, 0.05)) {
  s <- FindClusters(s, resolution = i, graph.name = 'harmony_snn')
}

# Plot clusters for each resolution
plt_list <- list()
index <- 1
for (i in seq(0.25, 1, 0.05)) {
  colname <- paste0("harmony_snn_res.", i)
  
  plt <- DimPlot(s, group.by = colname, reduction = "harmonyumap", raster = FALSE, shuffle = TRUE) + 
    scale_color_manual(values = randomcoloR::distinctColorPalette(length(unique(s@meta.data[,colname])))) +
    ggtitle(paste0("Resolution = ", i))
  
  plt_list[[index]] <- plt
  
  index <- index + 1
}
pdf(paste0(plot_dir, "combined_umap.pdf"), height = 60, width = 60)
ggpubr::ggarrange(plotlist = plt_list, nrow = 4, ncol = 4)
dev.off()

# Save plots individually
res <- seq(0.25, 1, 0.05)
for (index in 1:length(plt_list)) {
  pdf(paste0(plot_dir, "res_", res[[index]], ".pdf"), height = 15, width = 15)
  print(plt_list[[index]])
  dev.off()
}

# Save object
saveRDS(s, "/path/to/integrated/scRNAseq/object.rds")

