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
# Summary: Create Seurat objects using binned expression
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
})

# Define general output folder
output_folder <- "/path/to/preprocessing/output/folder/"

# Load python modules
reticulate::use_virtualenv("/path/to/r-miniconda/envs/r-reticulate")
scanpy <- reticulate::import('scanpy')
anndata <- reticulate::import('anndata')

# Initialize Seurat objects with binned data
for (sample in list.files(output_folder)) {
  adata <- scanpy$read_h5ad(paste0(output_folder, sample, "/binned_adata.h5ad"))
  adata$var_names_make_unique()
  counts <- adata$X
  colnames(counts) <- as.character(adata$var_names$values)
  rownames(counts) <- as.character(adata$obs$id)
  counts <- t(counts)
  s <- CreateSeuratObject(counts = counts)
  saveRDS(s, paste0(output_folder, sample, "/binned_seurat.rds"))
}

