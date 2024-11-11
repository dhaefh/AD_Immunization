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
# Summary: Preprocess RNA data with Seurat
#
#-------------------------------------------------------------------------------
# Initialization 

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

# Define paths to sample-level SpaceRanger output 
all_samples <- list.dirs(paste0(input_folder, "/data"), recursive = FALSE)

# Define paths to sample-level manual layer annotations
annotation_folder <- paste0(input_folder, "/layers/")

# Initialize list for sample objects 
objects <- list()

#-------------------------------------------------------------------------------
# Preprocess cohort 5/7 samples

# Initialize Seurat objects for each sample
i <- 1
for (cur_sample in all_samples) {
  print(cur_sample)
  cur_dir <- paste0(cur_sample, "/outs")
  cur_slice <- strsplit(cur_sample, "/")[[1]][[8]]
  
  # Replace dash with underscore
  cur_slice <- str_replace_all(cur_slice, "-", ".")
  
  # Load RNA data
  s <- Seurat::Load10X_Spatial(cur_dir, slice = cur_slice, assay = 'Spatial')
  
  # Convert coordinates to integer
  s@images[[cur_slice]]@coordinates[["tissue"]] <- as.integer(s@images[[cur_slice]]@coordinates[["tissue"]])
  s@images[[cur_slice]]@coordinates[["row"]] <- as.integer(s@images[[cur_slice]]@coordinates[["row"]])
  s@images[[cur_slice]]@coordinates[["col"]] <- as.integer(s@images[[cur_slice]]@coordinates[["col"]])
  s@images[[cur_slice]]@coordinates[["imagerow"]] <- as.integer(s@images[[cur_slice]]@coordinates[["imagerow"]])
  s@images[[cur_slice]]@coordinates[["imagecol"]] <- as.integer(s@images[[cur_slice]]@coordinates[["imagecol"]])
  
  # Log-normalize RNA count data
  s <- NormalizeData(s, assay = "Spatial")
  
  # Add spatial information to meta data
  cur_saptial_info <- paste0(cur_dir, "/spatial/tissue_positions_list.csv")
  print(cur_saptial_info)
  spatial_df <- read.csv(cur_saptial_info)
  row.names(spatial_df) <- spatial_df$barcode
  spatial_df <- spatial_df[match(row.names(s@meta.data), row.names(spatial_df)),]
  s <- AddMetaData(s, metadata = spatial_df)
  
  # Add manual annotations to meta data
  cur_anno_file <- paste0(annotation_folder, cur_slice, ".csv")
  anno_df <- read.csv(cur_anno_file)
  row.names(anno_df) <- anno_df$Barcode
  anno_df <- anno_df[match(row.names(s@meta.data), row.names(anno_df)),]
  s[['Manual_Layer']] <- anno_df$manual
  
  # Add sample ID to cell names
  s <- RenameCells(object = s, add.cell.id = cur_slice)
  s@meta.data$sample_id <- cur_slice
  s@meta.data$sample_barcode <- row.names(s@meta.data)
  
  objects[[i]] <- s
  i <- i + 1
}

# Merge all samples
all_samples_seurat <- merge(objects[[1]], objects[2:length(objects)], project = "cohort_5", merge.data = TRUE)

# Correct manual layer spelling
all_samples_seurat@meta.data$Manual_Layer[all_samples_seurat@meta.data$Manual_Layer == "non-hipppcampus"] <- "non-hippocampus"

# Save merged Seurat object
saveRDS(all_samples_seurat, glue("{output_folder}/all_samples_01_rna.rds"))