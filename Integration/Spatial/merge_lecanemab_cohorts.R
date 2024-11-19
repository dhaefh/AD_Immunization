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
# Summary: Merge all lecanemab data and perform final QC filtering based on manual layer and protein expression 
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
  library("stringr")
})

# Define output folder 
output_folder <- "/path/to/cohort578/integration/output/folder/"

# Define filter operator
`%notin%` <- Negate(`%in%`)

# Load post-QC cohort 5/7 objects
rna57 <- readRDS("/path/to/post/qc/cohort57/rna/object.rds")
pro57 <- readRDS("/path/to/post/qc/cohort57/protein/object.rds")

# Remove spots not located on cortical or hippocampal tissue
rna57 <- subset(rna57, layer_use %notin% c("non-hippocampus", "exclude"))
pro57 <- subset(pro57, cells = rownames(rna57@meta.data))

# Create temporary objects with merged sample data to identify spots with zero protein expression
temp_rna57 <- JoinLayers(rna57)
temp_pro57 <- JoinLayers(pro57)

# Remove spots with zero protein expression 
spots_remove <- colnames(temp_pro57[, colSums(temp_pro57@assays$Protein@layers$counts)==0])
rna57 <- subset(rna57, sample_barcode %notin% spots_remove)
pro57 <- subset(pro57, sample_barcode %notin% spots_remove)

# Save memory
temp_rna57 <- NULL
temo_pro57 <- NULL
gc()

# Load post-QC cohort 8 Seurat objects
rna8 <- readRDS("/path/to/post/qc/cohort8/rna/object.rds")
pro8 <- readRDS("/path/to/post/qc/cohort8/protein/object.rds")

# Remove spots not located on cortical or hippocampal tissue
rna8 <- subset(rna8, manual_layer %notin% c("non-hippocampus", ""))
pro8 <- subset(pro8, cells = rownames(rna8@meta.data))

# Create temporary objects with merged sample data to identify spots with zero protein expression
temp_rna8 <- JoinLayers(rna8)
temp_pro8  <- JoinLayers(pro8)

# Remove spots with zero protein expression 
spots_remove <- colnames(temp_pro8[, colSums(temp_pro8@assays$Protein@layers$counts)==0])
rna8 <- subset(rna8, sample_barcode %notin% spots_remove)
pro8 <- subset(pro8, sample_barcode %notin% spots_remove)

# Save memory
temp_rna8 <- NULL
temo_pro8 <- NULL
gc()

# Make sure meta data columns match 
rna57$Manual_Layer <- rna57$layer_adj
rna57@meta.data <- rna57@meta.data %>% rename(manual_layer = Manual_Layer)
rna57@meta.data <- rna57@meta.data[,1:15]
sum(colnames(rna57@meta.data) != colnames(rna8@meta.data))

# Merge cohorts and split layers by sample ID 
rna_merged <- merge(rna8, rna57)
pro_merged <- merge(pro8, pro57)
rna_merged <- JoinLayers(rna_merged)
rna_merged[['Spatial']] <- split(rna_merged[['Spatial']], f = rna_merged$sample_id)
pro_merged <- JoinLayers(pro_merged)
pro_merged[['Protein']] <- split(pro_merged[['Protein']], f = pro_merged$sample_id)

# Save objects (note that this is the final protein Seurat object used for downstream analysis)
saveRDS(rna_merged, paste0(output_folder, "rna_merged_for_integration.rds"))
saveRDS(pro_merged, paste0(output_folder, "pro_merged_for_integration.rds"))
 