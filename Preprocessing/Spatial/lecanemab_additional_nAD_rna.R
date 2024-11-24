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
# Summary: Preprocess Space Ranger output for additional lecanemab nAD samples with Seurat
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
input_folder <- "/path/to/spaceranger/output"
output_folder <- "/path/to/preprocessing/folder"

# Define paths to sample-level Space Ranger output 
all_samples <- list.dirs(input_folder, recursive = FALSE)

# Initialize list for sample objects
seurat_objects <- list()

# Initialize Seurat objects for each sample
i <- 1
for (cur_sample in all_samples) {
  
  # Define sample ID and path to load spatial data
  cur_dir <- paste0(cur_sample, "/outs")
  cur_slice <- strsplit(cur_sample, "/")[[1]][[length(strsplit(cur_sample, "/")[[1]])]]
  if (cur_slice == "A14-193_A1") {
    cur_slice <- str_replace(cur_slice, "_A1", "_1")
  }
  cur_slice <- str_replace(cur_slice, "Cohort8-B_", "")
  cur_slice <- str_replace(cur_slice, "Cohort8-A_", "")
  cur_slice <- str_replace_all(cur_slice, "-", ".")
  cur_slice <- str_replace_all(cur_slice, "_", ".")
  print(cur_slice)
  
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
  cur_saptial_info <- paste0(cur_dir, "/spatial/tissue_positions.csv")
  spatial_df <- read.csv(cur_saptial_info)
  row.names(spatial_df) <- spatial_df$barcode
  spatial_df <- spatial_df[row.names(cur_s@meta.data),]
  print(sum(rownames(spatial_df) != rownames(cur_s@meta.data)))
  print(sum(is.na(spatial_df)))
  cur_s <- AddMetaData(cur_s, metadata = spatial_df)
  
  # Add sample ID to cell names
  cur_s <- RenameCells(object = cur_s, add.cell.id = cur_slice)
  cur_s@meta.data$sample_id <- cur_slice
  cur_s@meta.data$sample_barcode <- row.names(cur_s@meta.data)
  
  # Update list of Seurat objects
  seurat_objects[[i]] <- cur_s
  i <- i + 1
}

# Merge sample objects
s <- merge(seurat_objects[[1]], seurat_objects[2:length(seurat_objects)], project = "Cohort_8", merge.data = TRUE)

# Save merged Seurat object
saveRDS(s, paste0(output_folder, "/s_preprocessed.rds"))
