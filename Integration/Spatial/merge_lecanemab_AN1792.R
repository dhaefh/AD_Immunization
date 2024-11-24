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
# Summary: Merge lecanemab and AN1792 data
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
  library("ggplot2")
})

# Define output folder 
output_folder <- "/path/to/all/cohorts/integration/output/folder/"

# Load integrated objects for AN1792 and lecanemab
cohort1 <- readRDS("/path/to/AN1792/integrated/object.rds")
cohort578 <- readRDS("/path/to/lecanemab/integrated/object.rds")

# Keep only raw counts in each object
DefaultAssay(cohort1) <- "Spatial"
DefaultAssay(cohort578) <- "Spatial"
cohort1 <- JoinLayers(cohort1)
cohort578 <- JoinLayers(cohort578)
cohort1 <- DietSeurat(cohort1, assays = "Spatial")
cohort578 <- DietSeurat(cohort578, assays = "Spatial")
gc()

# Make sure meta data columns match 
cohort578@meta.data[,c("nCount_Protein", "nFeature_Protein")] <- NULL
cohort1@meta.data <- cohort1@meta.data %>% rename(manual_layer = manual_annotation)
cohort578@meta.data[,c("cortical_amyloid", "vascular_amyloid")] <- NULL
cohort1@meta.data <- cohort1@meta.data[,colnames(cohort578@meta.data)]

# Merge objects
merged <- merge(cohort1, cohort578)
merged@meta.data[,c("nCount_SCT", "nFeature_SCT")] <- NULL
merged <- JoinLayers(merged)
merged[['Spatial']] <- split(merged[['Spatial']], f = merged$sample_id)
saveRDS(merged, paste0(output_folder, "all_merged_for_integration.rds"))

