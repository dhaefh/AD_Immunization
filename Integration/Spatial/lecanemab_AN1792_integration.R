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
suppressMessages ({
  library('Seurat')
  library('dplyr')
})

# Define output folder 
output_folder <- "/path/to/all/cohorts/integration/output/folder/"

# Load Seurat object 
s <- readRDS(paste0(output_folder, "all_merged_for_integration.rds"))

# Run SCTransform
s <- SCTransform(s, assay = 'Spatial')

# Run PCA
s <- RunPCA(s, verbose = F)

# Integrate data 
s <- IntegrateLayers(object = s, method = HarmonyIntegration, normalization.method = 'SCT', orig.reduction = 'pca', verbose = F, new.reduction = "harmony")

# Save object
saveRDS(s, paste0(output_folder, "all_merged_integrated.rds"))