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
# Summary: Clustering of lecanemab and AN1792 amyloid-rich spots in gray matter
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages ({
  library('Seurat')
  library('dplyr')
  library("tidyverse")
})

# Define output folder for plots
plot_dir <- "/path/to/cluster/plots/"

# Load all-cohorts gray matter amyloid-rich integrated Seurat object 
s <- readRDS("/path/to/integrated/all/cohorts/amyloid/object.rds")

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
  
  plt <- DimPlot(s, group.by = colname, reduction = "harmonyumap", raster = FALSE, shuffle = TRUE, pt.size = 0.75) + 
    scale_color_manual(values = randomcoloR::distinctColorPalette(length(unique(s@meta.data[,colname])))) +
    ggtitle(paste0("Resolution = ", i))
  
  plt_list[[index]] <- plt
  
  index <- index + 1
}
pdf(paste0(plot_dir, "combined_umap.pdf"), height = 40, width = 40)
ggpubr::ggarrange(plotlist = plt_list, nrow = 4, ncol = 4)
dev.off()

# Save plots individually
res <- seq(0.25, 1, 0.05)
for (index in 1:length(plt_list)) {
  pdf(paste0(plot_dir, "res_", res[[index]], ".pdf"), height = 10, width = 10)
  print(plt_list[[index]])
  dev.off()
}

# Save object
saveRDS(s, "/path/to/integrated/all/cohorts/amyloid/object.rds")
