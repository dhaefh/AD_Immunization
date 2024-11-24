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
# Summary: QC filtering for additional lecanemab nAD samples with Seurat
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
  library("ggplot2")
})

# Define filter operator
`%notin%` <- Negate(`%in%`)

# Define output folder 
output_folder <- "/path/to/preprocessing/folder"

# Define path to QC threshold file
threshold_file <- "/path/to/qc/metrics/file.xlsx"

# Load preprocessed Seurat object
s <- readRDS(paste0(output_folder, "/s_preprocessed.rds"))

# Load per-sample QC threshold data and set row names to sample ID
threshold_df <- readxl::read_xlsx(threshold_file) %>% data.frame()
colnames(threshold_df) <- c("sample_id", "max_umi", "min_umi", "min_features")
threshold_df$sample_id <- str_replace_all(threshold_df$sample_id, "-", ".")
threshold_df$sample_id <- str_replace_all(threshold_df$sample_id, "_", ".")
rownames(threshold_df) <- threshold_df$sample_id
sum(rownames(threshold_df) %notin% s$sample_id)

# Generate pre-QC spot count table
pre_qc_table <- as.data.frame(table(s$sample_id))
colnames(pre_qc_table)[1] <- "sample_id"
colnames(pre_qc_table)[2] <- "n_spots"
write.csv(pre_qc_table, paste0(output_folder, "/all_sample_pre_qc_spot_count.csv"), row.names = FALSE, quote = FALSE)

# Split into sample objects
s_list <- SplitObject(s, split.by = "sample_id")

# Initialize list for post-QC sample objects
objects_post_qc <- c()

# Initialize QC stats table
sample_qc_stats <- NULL

# QC filtering for each sample
for (cur_sample in names(s_list)) {
  
  # Get current sample object
  cur_s <- s_list[[cur_sample]]
  
  # Get current UMI/nFeature cutoffs (change column names if needed)
  cur_umi_min_cutoff <- as.numeric(threshold_df[cur_sample, "min_umi"])
  cur_umi_max_cutoff <- as.numeric(threshold_df[cur_sample, "max_umi"])
  cur_nfeature_min_cutoff <- as.numeric(threshold_df[cur_sample, "min_features"])
  
  # Pre-QC total spots for current sample
  cur_sample_spots_before_QC <- dim(cur_s)[2]
  
  # Remove spots based on UMI thresholds
  spots_to_discard_df <-cur_s@meta.data[which((cur_s$nCount_Spatial<cur_umi_min_cutoff) |
                                                (cur_s$nCount_Spatial>cur_umi_max_cutoff) |
                                                (cur_s$nFeature_Spatial<cur_nfeature_min_cutoff)),]
  spots_to_discard <- spots_to_discard_df|>rownames()
  cur_s <- subset(cur_s, subset = sample_barcode %notin% spots_to_discard)
  
  # Get number of spots removed based on UMI count
  cur_umi_spots <- length(spots_to_discard)
  
  # Remove edge spots
  max_row_idx <- max(cur_s$array_row)
  max_col_idx <- max(cur_s$array_col)
  min_row_idx <- min(cur_s$array_row)
  min_col_idx <- min(cur_s$array_col)
  spots_to_discard_df <- cur_s@meta.data[which((cur_s$array_col==min_col_idx) | (cur_s$array_col==max_col_idx) | 
                                                 (cur_s$array_row==max_row_idx) | (cur_s$array_row==min_row_idx)),]
  spots_to_discard <- spots_to_discard_df|>rownames()
  cur_s <- subset(cur_s, subset = sample_barcode %notin% spots_to_discard)
  
  # Get number of edge spots removed
  cur_edge_spots <- length(spots_to_discard)
  
  # Filter using MT % thresholds (30% for HIPP, 20% for all other regions)
  cur_s[["percent.mt"]] <- PercentageFeatureSet(cur_s, pattern = "^MT-", assay = 'Spatial')
  pre_MT <- nrow(cur_s@meta.data)
  if (cur_sample %in% c("A14.193.9", "A11.170.9")) {
    print(paste0("Filtering ", cur_sample, " using 30% threshold"))
    cur_s <- subset(cur_s, subset = percent.mt < 30)
  } else {
    cur_s <- subset(cur_s, subset = percent.mt < 20)
  }
  post_MT <- nrow(cur_s@meta.data)
  
  # Get number of spots removed based on MT %
  cur_mt_spots <- pre_MT - post_MT
  
  # Update object list
  objects_post_qc <- c(objects_post_qc, cur_s)
  
  # Form QC tables
  cur_sample_spots_after_QC <- dim(cur_s)[2]
  cur_sample_stat_df <- data.frame(list(
    "sample_id" = cur_sample,
    "spots_before_QC" = cur_sample_spots_before_QC,
    "spots_after_QC" = cur_sample_spots_after_QC,
    "total_spots_removed" = cur_sample_spots_before_QC - cur_sample_spots_after_QC,
    "umi_nfeat_outliers" = cur_umi_spots,
    "percent_umi_nfeat_outliers" = (cur_umi_spots/cur_sample_spots_before_QC)*100,
    "edge_spots" = cur_edge_spots,
    "percent_edge_spots" = (cur_edge_spots/cur_sample_spots_before_QC)*100,
    "MT_spots" = cur_mt_spots,
    "percent_MT_spots" = (cur_mt_spots/cur_sample_spots_before_QC)*100
  ))
  sample_qc_stats <- rbindlist(list(
    sample_qc_stats, cur_sample_stat_df
  ))|>
    data.frame()
}

# Save per-sample QC stats
write.csv(sample_qc_stats, paste0(output_folder, "/all_sample_qc_stats.csv"), row.names = FALSE, quote = FALSE)

# Merge sample objects 
s <- merge(objects_post_qc[[1]], objects_post_qc[2:length(objects_post_qc)], merge.data = TRUE)

# Split layers by sample ID 
s <- JoinLayers(s)
s[['Spatial']] <- split(s[['Spatial']], f = s$sample_id)

# Save object
saveRDS(s, file = paste0(output_folder, "/s_post_qc.rds"))

# Filter protein object to post-QC spots
pro <- readRDS("/path/to/preprocessed/protein/object.rds")
pro <- subset(pro, cells = rownames(s@meta.data))
saveRDS(pro, "/path/to/post/qc/protein/object.rds")
