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
# Summary: Preprocess Space Ranger output for lecanemab with Seurat
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library('Seurat')
  library('glue')
  library('dplyr')
  library('stringr')
  library("Matrix")
  library("hdf5r")
})

# Define input and output folders
input_folder <- "/path/to/input/data"
output_folder <- "/path/to/output/folder/"

# Define paths to sample-level Space Ranger output 
all_samples <- list.dirs(paste0(input_folder, "/data"), recursive = FALSE)

# Initialize list for sample objects 
objects <- list()

# Initialize Seurat objects for each sample
i <- 1
for (cur_sample in all_samples) {
  print(cur_sample)
  cur_dir <- paste0(cur_sample, "/outs")
  cur_slice <- strsplit(cur_sample, "/")[[1]][[8]]
  
  # Format sample ID
  cur_slice <- str_replace_all(cur_slice, "-", ".")
  
  # Load protein data (isotype-normalized raw counts)
  s <- Seurat::Load10X_Spatial(cur_dir, slice = cur_slice, assay = 'Protein', filename = "filtered_feature_bc_matrix.h5")
  
  # Remove isotype controls
  features <- rownames(s)[-(grep("^mouse|^rat", rownames(s)))]
  
  # Ensure that protein assay is v5
  pro_assay <- CreateAssay5Object(counts = s@assays$Protein@counts[features,])
  s[['Protein']] <- pro_assay
  
  # Convert coordinates to integer
  s@images[[cur_slice]]@coordinates[["tissue"]] <- as.integer(s@images[[cur_slice]]@coordinates[["tissue"]])
  s@images[[cur_slice]]@coordinates[["row"]] <- as.integer(s@images[[cur_slice]]@coordinates[["row"]])
  s@images[[cur_slice]]@coordinates[["col"]] <- as.integer(s@images[[cur_slice]]@coordinates[["col"]])
  s@images[[cur_slice]]@coordinates[["imagerow"]] <- as.integer(s@images[[cur_slice]]@coordinates[["imagerow"]])
  s@images[[cur_slice]]@coordinates[["imagecol"]] <- as.integer(s@images[[cur_slice]]@coordinates[["imagecol"]])
  
  # Add sample ID to cell names
  s <- RenameCells(object = s, add.cell.id = cur_slice)
  s@meta.data$sample_id <- cur_slice
  s@meta.data$sample_barcode <- row.names(s@meta.data)
  
  objects[[i]] <- s
  i <- i + 1
}

# Merge sample objects 
all_samples_seurat <- merge(objects[[1]], objects[2:length(objects)], project = "cohort_5", merge.data = TRUE)

# Save merged Seurat object
saveRDS(all_samples_seurat, paste0(output_folder, "all_samples_01_pro_filtered.rds"))