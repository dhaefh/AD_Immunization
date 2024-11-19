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
# Summary: Define microglia cluster enrichment using Cell2Location predicted abundance
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
  library("tidyr")
})

# Load C2L q05 cell type abundance predictions
meta <- read.csv("/path/to/cell2location/mapping/q05_meta.csv", row.names = 1)

# Define number of standard deviations above the mean to use for enrichment threshold and create output folder
n_sd <- 3
output_dir <- paste0("/path/to/cell2location/mapping/mg_enrichment/", n_sd, "_sd/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Calculate mean and standard deviation per microglia cluster per sample 
summary <- meta[,c("sample_id", paste0("Mg.", 0:4))]
summary <- summary %>% pivot_longer(cols = paste0("Mg.", 0:4), names_to = "cluster", values_to = "q05")
summary <- summary %>% group_by(sample_id, cluster) %>% summarize(mean = mean(q05), sd = sd(q05))
summary$enrichment_thresh <- summary$mean + n_sd*summary$sd

# Save summary data
write.csv(summary, paste0(output_dir, "sample_thresholds.csv"), row.names = FALSE)

# Define enrichment for each microglia cluster using sample specific threshold of n standard deviations above the mean 
for (clust in paste0("Mg.", 0:4)) {
  
  # Initialize enrichment variable
  meta[[paste0(clust, "_enriched")]] <- 0
  
  # Define enrichment in each sample
  for (sample in unique(meta$sample_id)) {
    thresh <- summary$enrichment_thresh[summary$sample_id == sample & summary$cluster == clust]
    meta[[paste0(clust, "_enriched")]][meta[[clust]] > thresh & meta$sample_id == sample] <- 1
  }
}

# Save meta data with enrichment labels
write.csv(meta, paste0(output_dir, "q05_meta_with_enrichment.csv"))






