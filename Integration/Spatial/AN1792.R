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
# Summary: Run SCTransform and perform CCA integration in Seurat
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages ({
  library('Seurat')
  library('glue')
  library('dplyr')
  library('stringr')
  library("Matrix")
  library("data.table")
  library('sctransform')
  library('ggplot2')
  library('glmGamPoi')
})

# Define input/output folders
prev_out_dir <- "/path/to/qc/output/folder"
objects_out_dir <- "/path/to/data/output/folder" 

# Load post-QC Seurat object
all_seurat_s <- glue("{prev_out_dir}/all_samples_02.rds")|>readRDS()

# Remove spots not covering tissue
all_seurat_s <- subset(all_seurat_s, manual_annotation != "")

# Run SCTransform
all_seurat_s <- SCTransform(all_seurat_s, assay = 'Spatial')

# Run PCA
all_seurat_s <- RunPCA(all_seurat_s, npcs = 30, verbose = F)

# Integrate data 
all_seurat_s <- IntegrateLayers(object = all_seurat_s, method = CCAIntegration, normalization.method = 'SCT', verbose = F)

# Save object 
saveRDS(all_seurat_s, paste0(objects_out_dir, "/all_samples_03.rds"))