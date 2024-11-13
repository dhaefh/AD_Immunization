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
# Summary: QC filtering with Seurat
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages ({
  library('Seurat')
  library('glue')
  library('dplyr')
  library('stringr')
  library("Matrix")
  library("data.table")
  library("ie2misc")
  library("ggplot2")
  library('SeuratDisk')
})

# Define filter operator
`%notin%` <- Negate(`%in%`)

# Define input/output folders and filepaths
prev_out_dir <- "/path/to/preprocessing/folder"
all_seurat_file <- glue("{prev_out_dir}/all_samples_01_rna.rds")
objects_out_dir <- "/path/to/objects/output/folder"
QC_per_sample_metric_file <- "/path/to/qc/metrics/file.xlsx"

#-------------------------------------------------------------------------------
# QC filtering for Lecanemab cohort 5/7 

# Read per sample QC metric file
QC_per_sample_metric_df <- readxl::read_xlsx(QC_per_sample_metric_file)|>as.data.frame()
QC_per_sample_metric_df[["sample_id"]] <- sapply(QC_per_sample_metric_df[["file-name"]], function(x){str_replace_all(x,"-", ".")})
rownames(QC_per_sample_metric_df) <- QC_per_sample_metric_df[["sample_id"]]

# Initialize variable for per-sample QC stats
sample_qc_stats <- NULL

# Load preprocessed Seurat object
all_seurat <- readRDS(all_seurat_file)

# Generate pre-QC spot count table
pre_qc_table <- all_seurat$sample_id|>table()|>as.data.frame()
colnames(pre_qc_table)[1] <- "sample_id"
colnames(pre_qc_table)[2] <- "n_spots"
QC_stats_out_file <- glue("{objects_out_dir}/all_sample_pre_qc_spot_count.csv")
write.csv(pre_qc_table, QC_stats_out_file, row.names = FALSE, quote = FALSE)

# Split merged Seurat object into sample objects, initialize list for post-QC objects
s_list <- SplitObject(all_seurat, split.by = "sample_id")
post_qc <- c()

# Apply QC filtering to each sample object
for (cur_sample in names(s_list)) {
  
  cur_s <- s_list[[cur_sample]]
  
  # Get UMI and nFeature thresholds for current sample
  cur_umi_min_cutoff <- QC_per_sample_metric_df[cur_sample, "min-UMIs-per-barcode"]
  cur_umi_max_cutoff <- QC_per_sample_metric_df[cur_sample, "max-UMIs-per-barcode"]
  cur_nfeature_min_cutoff <- QC_per_sample_metric_df[cur_sample, "min-genes-per-barcode"]
  
  # Remove images from other samples
  all_images <- cur_s@images|>names()
  images_to_remove <- all_images[-grep(cur_sample, all_images, fixed=T)]
  for(cur_img_to_remove in images_to_remove){
    cur_s@images[[cur_img_to_remove]] <- NULL
  }
  
  # Pre-QC total spots for current sample
  cur_sample_spots_before_QC <- dim(cur_s)[2]
  
  # Remove spots based on UMI/nFeature thresholds
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
  
  # Remove spots based on MT %
  cur_s[["percent.mt"]] <- PercentageFeatureSet(cur_s, pattern = "^MT-", assay = 'Spatial')
  pre_MT <- nrow(cur_s@meta.data)
  
  # Sample-specific MT thresholds
  if (unique(cur_s@meta.data$sample_id %in% c("NMA22.A9", "NMA22.B9"))) {
    cur_s <- subset(cur_s, subset = percent.mt < 30)
  } else {
    cur_s <- subset(cur_s, subset = percent.mt < 20)
  }
  
  post_MT <- nrow(cur_s@meta.data)
  
  # Get number of spots removed based on MT %
  cur_mt_spots <- pre_MT - post_MT
  
  post_qc <- c(post_qc, cur_s)
  
  # Update QC stats table
  cur_sample_spots_after_QC <- dim(cur_s)[2]
  cur_sample_stat_df <- data.frame(list(
    "sample_id" = cur_sample,
    "spots_before_QC" = cur_sample_spots_before_QC,
    "spots_after_QC" = cur_sample_spots_after_QC,
    "total_spots_removed" = cur_sample_spots_before_QC - cur_sample_spots_after_QC,
    "umi_outliers" = cur_umi_spots,
    "percent_umi_outliers" = (cur_umi_spots/cur_sample_spots_before_QC)*100,
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

# Save table of QC stats
sample_QC_stats_out_file <- glue("{objects_out_dir}/all_sample_post_qc_spot_count.csv")
write.csv(sample_qc_stats, sample_QC_stats_out_file, row.names = FALSE, quote = FALSE)

# Merge post-QC sample objects
all_seurat <- merge(post_qc[[1]], post_qc[2:length(post_qc)], merge.data = TRUE)

# Filter protein objects to post-QC spots
pro_filtered <- readRDS(glue("{prev_out_dir}/all_samples_01_pro_filtered.rds"))
pro_filtered <- subset(pro_filtered, sample_barcode %in% all_seurat@meta.data$sample_barcode)

# Save objects
saveRDS(all_seurat, file = glue("{objects_out_dir}/all_samples_02_rna.rds"))
saveRDS(pro_filtered, file = glue("{objects_out_dir}/all_samples_02_pro_filtered.rds"))