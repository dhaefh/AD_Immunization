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
# Summary: Generate LOESS predictions for iAD gene expression based on amyloid density in the plaque niche
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ggrepel")
  library("ggthemes")
  library("ggpubr")
  library("Seurat")
  library("factoextra")
  library("scales")
  library("patchwork")
  library("clusterProfiler")
  library("pheatmap")
  library("colorspace")
  library("stringi")
  library("randomcoloR")
})

# Define output folder
output_folder <- "/path/to/AN1792/loess/data/"

# Define filter operator 
`%notin%` <- Negate(`%in%`)

# Load integrated AN1792 Seurat object
s <- readRDS("/path/to/integrated/AN1792/object.rds")

# Subset for iAD and nAD cortical amyloid-rich spots and first + second order neighbors in gray matter
s <- subset(s, amyloid_neighbor_final %in% c("amyloid", "first_neighbor", "second_neighbor") & manual_annotation != "white" & condition %in% c("iAD", "nAD"))
gc()

# Recorrect SCT data
for (name in names(s@assays$SCT@SCTModel.list)) {
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}
s <- PrepSCTFindMarkers(object = s)

# Extract recorrected SCT expression and meta data from Seurat object
all_data <- GetAssayData(s, assay = "SCT", layer = "data")
meta <- s@meta.data

# Save memory 
s <- NULL
gc()

# Convert expression matrix to data frame
data <- data.frame(all_data)

# Save memory 
all_data <- NULL
gc()

# Format data for LOESS
data <- t(data)
data <- data.frame(data) 
data$row_name <- stri_replace_last_fixed(row.names(data), ".", "-")
row.names(data) <- data$row_name
print(sum(row.names(data) != row.names(meta)))
data$amyloid <- meta$amyloid_fluo
data$condition <- meta$condition

# Define function to generate LOESS predictions 
run_loess <- function(gene, data) {
  
  # Extract gene expression and amyloid density from input data
  data <- data[, c(gene, "amyloid")]
  
  # Fit LOESS model with amyloid density as predictor for gene expression
  lo <- loess(get(gene) ~ amyloid, data)
  
  # Generate LOESS predictions for a uniform grid of amyloid density values
  lo_predict <- predict(lo, data.frame(amyloid = seq(min(data$amyloid), max(data$amyloid), 1))) %>% as.data.frame()
  colnames(lo_predict) <- gene
  
  return(lo_predict)
}

# Ensure names of genes to test exactly match column names 
gene_names <- colnames(data)[colnames(data) %notin% c("row_name", "amyloid", "condition")]

# Generate predictions for iAD
cur_data <- data[data$condition == "iAD", c(gene_names, "amyloid")]
predictions <- lapply(gene_names, run_loess, data = cur_data)
merged <- as.data.frame(predictions)
merged$amyloid <- seq(min(cur_data$amyloid), max(cur_data$amyloid), 1)
merged <- dplyr::relocate(merged, amyloid)
write.csv(merged, paste0(output_folder, "iAD_all_gene_predictions_broad.csv"), row.names = FALSE)

