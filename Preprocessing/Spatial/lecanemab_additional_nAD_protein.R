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
objects <- list()

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
  
  # Load protein data (isotype-normalized raw counts)
  cur_s <- Seurat::Load10X_Spatial(cur_dir, slice = cur_slice, assay = 'Protein', filename = "filtered_feature_bc_matrix.h5")
  
  # Remove isotype controls
  features <- rownames(cur_s)[-(grep("^mouse|^rat", rownames(cur_s)))]
  
  # Ensure that protein assay is v5
  pro_assay <- CreateAssay5Object(counts = cur_s@assays$Protein@counts[features,])
  cur_s[['Protein']] <- pro_assay
  
  # Convert coordinates to integer
  cur_s@images[[cur_slice]]@coordinates[["tissue"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["tissue"]])
  cur_s@images[[cur_slice]]@coordinates[["row"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["row"]])
  cur_s@images[[cur_slice]]@coordinates[["col"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["col"]])
  cur_s@images[[cur_slice]]@coordinates[["imagerow"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["imagerow"]])
  cur_s@images[[cur_slice]]@coordinates[["imagecol"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["imagecol"]])
  
  cur_s <- RenameCells(object = cur_s, add.cell.id = cur_slice)
  cur_s@meta.data$sample_id <- cur_slice
  cur_s@meta.data$sample_barcode <- row.names(cur_s@meta.data)
  
  objects[[i]] <- cur_s
  i <- i + 1
}

# Merge sample objects 
s <- merge(objects[[1]], objects[2:length(objects)], project = "Cohort_8")

# Save merged Seurat object
saveRDS(s, paste0(output_folder, "/s_preprocessed.rds"))



