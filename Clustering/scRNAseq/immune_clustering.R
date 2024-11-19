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
# Summary: Clustering of scRNAseq immune cells in Seurat
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
})

# Define output folder
output_dir <- "/path/to/clustering/output/"

# Create subfolders for plots
plot_dir <- paste0(output_dir, "plots/cca_immune_pool_merged/")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(plot_dir, "20_pcs"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(plot_dir, "30_pcs"), showWarnings = FALSE, recursive = TRUE)

# Load integrated scRNAseq immune Seurat object
s <- readRDS("/path/to/integrated/scRNAseq/immune/object.rds")

# Find clusters for a sequence of resolutions, using 20 and 30 PCs (30 PCs are ultimately used to define clusters)
s <- FindNeighbors(s, reduction = "cca", dims = 1:20, graph.name = c('cca_nn20pcs', 'cca_snn20pcs'))
for (i in seq(0.05, 1, 0.05)) {
  s <- FindClusters(s, resolution = i, graph.name = 'cca_snn20pcs')
}
s <- FindNeighbors(s, reduction = "cca", dims = 1:30, graph.name = c('cca_nn30pcs', 'cca_snn30pcs'))
for (i in seq(0.05, 1, 0.05)) {
  s <- FindClusters(s, resolution = i, graph.name = 'cca_snn30pcs')
}

# Save object 
saveRDS(s, "/path/to/integrated/scRNAseq/immune/object.rds")

# Plot clusters using 20 PCs
plt_list <- list()
index <- 1
for (i in seq(0.05, 1, 0.05)) {
  colname <- paste0("cca_snn20pcs_res.", i)
  
  plt <- DimPlot(s, group.by = colname, reduction = "ccaumap20pcs", raster = FALSE, shuffle = TRUE, pt.size = 2) + 
    scale_color_manual(values = randomcoloR::distinctColorPalette(length(unique(s@meta.data[,colname])))) +
    ggtitle(paste0("Resolution = ", i))
  
  plt_list[[index]] <- plt
  
  index <- index + 1
}
pdf(paste0(plot_dir, "20_pcs/combined_umap.pdf"), height = 60, width = 75)
ggpubr::ggarrange(plotlist = plt_list, nrow = 4, ncol = 5)
dev.off()

# Save plots individually
res <- seq(0.05, 1, 0.05)
for (index in 1:length(plt_list)) {
  pdf(paste0(plot_dir, "20_pcs/res_", res[[index]], ".pdf"), height = 15, width = 15)
  print(plt_list[[index]])
  dev.off()
}

# Plot clusters using 30 PCs
plt_list <- list()
index <- 1
for (i in seq(0.05, 1, 0.05)) {
  colname <- paste0("cca_snn30pcs_res.", i)
  
  plt <- DimPlot(s, group.by = colname, reduction = "ccaumap30pcs", raster = FALSE, shuffle = TRUE, pt.size = 2) + 
    scale_color_manual(values = randomcoloR::distinctColorPalette(length(unique(s@meta.data[,colname])))) +
    ggtitle(paste0("Resolution = ", i))
  
  plt_list[[index]] <- plt
  
  index <- index + 1
}
pdf(paste0(plot_dir, "30_pcs/combined_umap.pdf"), height = 60, width = 75)
ggpubr::ggarrange(plotlist = plt_list, nrow = 4, ncol = 5)
dev.off()

# Save plots individually
res <- seq(0.05, 1, 0.05)
for (index in 1:length(plt_list)) {
  pdf(paste0(plot_dir, "30_pcs/res_", res[[index]], ".pdf"), height = 15, width = 15)
  print(plt_list[[index]])
  dev.off()
}

