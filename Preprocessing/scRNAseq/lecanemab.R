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
# Summary: Remove background contamination for lecanemab with SoupX 
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("Matrix")
  library("SoupX")
  library("DropletUtils")
  library("purrr")
  library("Hmisc")
})

# Define general output folder
general_dir <- "/path/to/lecanemab/SoupX/output/"

# Define function to run SoupX
run_soupx <- function(dir) {
  sample <- unlist(strsplit(dir, "/")) %>%
    tail(1)
  
  sample_dir <- paste0(data_dir, sample, "/")
  dir.create(sample_dir, showWarnings = FALSE)
  
  rounded_dir <- paste0(sample_dir, "rounded/")
  dir.create(rounded_dir, showWarnings = FALSE)
  
  unrounded_dir <- paste0(sample_dir, "unrounded/")
  dir.create(unrounded_dir, showWarnings = FALSE)
  
  print(paste0("Processing sample: ", sample))
  
  # Initialize Seurat object
  seurat <- Read10X(paste0(dir, "/count/sample_filtered_feature_bc_matrix")) %>%
    CreateSeuratObject()
  
  # Run SCTransform and PCA 
  seurat <- SCTransform(seurat) %>% RunPCA()
  
  # Generate clusters using default parameters (SoupX is not sensitive to clustering quality - Seurat default should be fine)
  seurat <- seurat %>% FindNeighbors(dims = 1:10) %>% FindClusters() 
  
  # Load in 10X data for SoupX 
  toc = Seurat::Read10X(paste0(dir,"/count/sample_filtered_feature_bc_matrix"))
  tod = Seurat::Read10X(paste0(dir,"/count/sample_raw_feature_bc_matrix"))
  
  # Filter table of droplets (tod) so that it has same number of genes as table of counts (toc)
  tod <- tod[which(rownames(tod) %in% rownames(toc)),]
  
  # Create SoupX object (default calcSoupProfile = TRUE)
  soupx = SoupChannel(tod, toc)
  
  # Add cluster info to SoupX object
  soupx <- setClusters(soupx, setNames(seurat$seurat_clusters, colnames(seurat)))
  
  # Save automated estimate plot of contamination fraction 
  pdf(file = paste0(plot_dir, sample, "_rho_distribution.pdf"), width = 4, height = 4)
  soupx <- autoEstCont(soupx, forceAccept = TRUE) # Set forceAccept = TRUE to allow to proceed with high contamination fraction, can look into problematic samples later and decide to exclude or manually set the fraction more in line with other samples
  dev.off()
  
  # Adjust counts, rounded and unrounded to integer
  print(paste0(sample, " Contamination Fraction: ", soupx$metaData$rho[1])) # Print contamination fraction
  rounded <- adjustCounts(soupx, roundToInt = TRUE)
  unrounded <- adjustCounts(soupx)

  # Save corrected counts
  DropletUtils:::write10xCounts(rounded_dir, rounded, overwrite = TRUE)
  DropletUtils:::write10xCounts(unrounded_dir, unrounded, overwrite = TRUE)
  
  return(soupx$metaData$rho[1])
}

# Define function to get sample names corresponding to contamination fractions
get_samples <- function(dir) {
  sample <- unlist(strsplit(dir, "/")) %>% tail(1)
  return(sample)
}

# Define paths to Cell Ranger output for each pool
cellranger_dirs <- list.dirs("/path/to/cellranger/output", full.names = TRUE, recursive = FALSE) 

# Assign pool names
names(cellranger_dirs) <- c("CAA", "pool1", "pool2")

# Run SoupX for each sample in each pool
for (pool in names(cellranger_dirs)) {
  
  # Get path to cellranger output for current pool 
  cur_cellranger_dir <- cellranger_dirs[[pool]]
  cur_cellranger_dir <- paste0(cur_cellranger_dir, "/outs/per_sample_outs")
  
  # Create general output folder
  output_dir <- paste0(general_dir, pool, "/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create plot subfolder
  plot_dir <- paste0(output_dir, "plots/")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create subfolder for corrected counts
  data_dir <- paste0(output_dir, "data/")
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  
  # List paths to sample-level cellranger output 
  sample_dirs <- list.dirs(cur_cellranger_dir, recursive = FALSE)
  print(sample_dirs)
  
  # Run SoupX
  rho <- lapply(sample_dirs, run_soupx)
  
  # Export contamination fractions per sample
  samples <- lapply(sample_dirs, get_samples)
  summary <- data.frame(sample = unlist(samples), rho = unlist(rho))
  write.csv(summary, file = paste0(output_dir, "/contamination_fractions.csv"))
}








