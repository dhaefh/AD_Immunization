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
# Summary: Preprocess Space Ranger output for AN1792 with Seurat
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
output_folder <- "/path/to/preprocessing/folder"

# Initialize list for sample objects
seurat_objects <- list()

# Define paths to sample-level Space Ranger output for original cohort
all_samples <- list.dirs(paste0(input_folder, "/cohort_1_data"), recursive = FALSE)

# Initialize Seurat objects for each sample
i <- 1
for (cur_sample in all_samples) {
  print(cur_sample)
  cur_dir <- paste0(cur_sample, "/outs")
  cur_slice_raw <- strsplit(cur_sample, "/")[[1]][[7]]
  
  # Replace dash with underscore
  cur_slice <- str_replace_all(cur_slice_raw, "-", ".")
  
  # Initialize Seurat object
  cur_s <- Seurat::Load10X_Spatial(cur_dir, slice = cur_slice, assay = 'Spatial')
  
  # Convert coordinates to integer
  cur_s@images[[cur_slice]]@coordinates[["tissue"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["tissue"]])
  cur_s@images[[cur_slice]]@coordinates[["row"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["row"]])
  cur_s@images[[cur_slice]]@coordinates[["col"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["col"]])
  cur_s@images[[cur_slice]]@coordinates[["imagerow"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["imagerow"]])
  cur_s@images[[cur_slice]]@coordinates[["imagecol"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["imagecol"]])
  
  # Log-normalize counts
  cur_s <- NormalizeData(cur_s, assay = "Spatial")
  
  # Add spatial information to meta data
  cur_saptial_info <- paste0(cur_dir, "/spatial/tissue_positions_list.csv")
  print(cur_saptial_info)
  spatial_df <- read.csv(cur_saptial_info)
  row.names(spatial_df) <- spatial_df$barcode
  spatial_df <- spatial_df[match(row.names(cur_s@meta.data), row.names(spatial_df)),]
  cur_s <- AddMetaData(cur_s, metadata = spatial_df)
  
  # Add sample ID to cell names
  cur_s <- RenameCells(object = cur_s, add.cell.id = cur_slice)
  cur_s@meta.data$sample_id <- cur_slice
  cur_s@meta.data$sample_barcode <- row.names(cur_s@meta.data)
  
  seurat_objects[[i]] <- cur_s
  i <- i + 1
  
}

# Define paths to sample-level Space Ranger output for additional NNC samples
all_samples <- list.dirs(paste0(input_folder, "/cohort_6_data"), recursive = FALSE)

# Initialize Seurat objects for each sample
for (cur_sample in all_samples) {
  print(cur_sample)
  cur_dir <- paste0(cur_sample, "/outs")
  cur_slice_raw <- strsplit(cur_sample, "/")[[1]][[7]]
  
  # Replace dash with underscore
  cur_slice <- str_replace_all(cur_slice_raw, "-", ".")
  
  # Initialize Seurat object
  cur_s <- Seurat::Load10X_Spatial(cur_dir, slice = cur_slice, assay = 'Spatial')
  
  # Convert coordinates to integer
  cur_s@images[[cur_slice]]@coordinates[["tissue"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["tissue"]])
  cur_s@images[[cur_slice]]@coordinates[["row"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["row"]])
  cur_s@images[[cur_slice]]@coordinates[["col"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["col"]])
  cur_s@images[[cur_slice]]@coordinates[["imagerow"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["imagerow"]])
  cur_s@images[[cur_slice]]@coordinates[["imagecol"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["imagecol"]])
  
  # Log-normalize counts
  cur_s <- NormalizeData(cur_s, assay = "Spatial")
  
  # Add spatial information to meta data
  cur_saptial_info <- paste0(cur_dir, "/spatial/tissue_positions_list.csv")
  print(cur_saptial_info)
  spatial_df <- read.csv(cur_saptial_info)
  row.names(spatial_df) <- spatial_df$barcode
  spatial_df <- spatial_df[match(row.names(cur_s@meta.data), row.names(spatial_df)),]
  cur_s <- AddMetaData(cur_s, metadata = spatial_df)
  
  # Add sample ID to cell names
  cur_s <- RenameCells(object = cur_s, add.cell.id = cur_slice)
  cur_s@meta.data$sample_id <- cur_slice
  cur_s@meta.data$sample_barcode <- row.names(cur_s@meta.data)
  
  seurat_objects[[i]] <- cur_s
  i <- i + 1
}

# Merge all samples
all_samples_seurat <- merge(seurat_objects[[1]], seurat_objects[2:length(seurat_objects)], project = "cohort_1_all", merge.data = TRUE)

# Save merged Seurat object
saveRDS(all_samples_seurat, glue("{output_folder}/all_samples_01.rds"))