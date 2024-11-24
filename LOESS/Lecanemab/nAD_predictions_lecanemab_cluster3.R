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
# Summary: Generate LOESS predictions for nAD gene expression based on amyloid density in the plaque niche, for genes in lecanemab cluster 3
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
output_folder <- "/path/to/lecanemab/loess/data/"

# Define filter operator 
`%notin%` <- Negate(`%in%`)

# Load LCMB cluster 3 genes
clusters <- read.csv(paste0(output_folder, "LCMB_clusters_all_genes_broad.csv"), row.names = 1)
clusters$gene <- row.names(clusters)
clust_genes <- clusters$gene[clusters$nclust_11 == 3]
clust_genes[grep("^KRTAP|^MT|^ERV", clust_genes)] <- str_replace(clust_genes[grep("^KRTAP|^MT|^ERV", clust_genes)], "\\.", "-")

# Load integrated lecanemab Seurat object
s <- readRDS("/path/to/integrated/lecanemab/object.rds")

# Filter for amyloid-rich spots and first + second order neighbors in gray matter
gray_layers <- unique(s$manual_layer[grep("gray", s$manual_layer)])
s <- subset(s, amyloid_neighbor_final %in% c("amyloid", "first_neighbor", "second_neighbor") & manual_layer %in% gray_layers)
gc()

# Recorrect SCT data 
for (name in names(s@assays$SCT@SCTModel.list)) {
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}
s <- PrepSCTFindMarkers(object = s)

# Extract recorrected SCT expression and meta data from Seurat object
all_data <- GetAssayData(s, assay = "SCT", layer = "data")
all_data <- all_data[row.names(all_data) %in% clust_genes,]
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
data$amyloid <- meta$cortical_amyloid
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
clust_genes <- colnames(data)[colnames(data) %notin% c("row_name", "amyloid", "condition")]
print(length(clust_genes))

# Generate predictions for CAA
cur_data <- data[data$condition == "CAA", c(clust_genes, "amyloid")]
predictions <- lapply(clust_genes, run_loess, data = cur_data)
merged <- as.data.frame(predictions)
merged$amyloid <- seq(min(cur_data$amyloid), max(cur_data$amyloid), 1)
merged <- dplyr::relocate(merged, amyloid)
write.csv(merged, paste0(output_folder, "CAA_predictions_clust3.csv"), row.names = FALSE)

