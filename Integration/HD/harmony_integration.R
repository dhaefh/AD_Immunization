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
# Summary: Run SCTransform and perform Harmony integration in Seurat
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
output_folder <- "/path/to/integration/output/folder/"

# Load merged post-QC Seurat object
s <- readRDS("/path/to/qc/output/folder/data/s_post_qc.rds")

# Run SCTransform 
s <- SCTransform(s)

# Run PCA
s <- RunPCA(s, verbose = F)

# Integrate data 
s <- IntegrateLayers(object = s, method = HarmonyIntegration, normalization.method = 'SCT', orig.reduction = 'pca', verbose = F, new.reduction = "harmony")

# Generate UMAP
ElbowPlot(s, ndims = 50)
s <- RunUMAP(s, dims = 1:30, reduction.name = 'harmonyumap', reduction = 'harmony')

# Save object
saveRDS(s, paste0(output_folder, "data/s_integrated.rds"))

