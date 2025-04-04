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
# Summary: Pre-DoubletFinder QC filtering for AN1792 with Seurat
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("ggthemes")
  library("ggrepel")
  library("grid")
  library("DoubletFinder")
  library("doMC")
  library("xlsx")
  library("RColorBrewer")
})

# Define output folder
output_dir <- "/path/to/AN1792/QC/output/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(output_dir, "data"), recursive = TRUE, showWarnings = FALSE)

# Define function to load Seurat objects (set rounded = FALSE to load non-integer SoupX counts)
load_seurat <- function(dir, id, rounded = TRUE) {
  
  if (rounded) {
    counts <- Read10X(paste0(dir, "/rounded"))
  } else {
    counts <- Read10X(paste0(dir, "/unrounded"))
  }
  
  return(CreateSeuratObject(counts = counts, project = id)) 
}

# Initialize Seurat objects for each pool using SoupX corrected counts
i <- 1
pool_list <- list()
for (pool in list.files("/path/to/AN1792/SoupX/output")) {
  soupx_dir <- paste0("/path/to/AN1792/SoupX/output/", pool, "/data")
  
  # List paths to SoupX corrected counts for each sample
  soupx_sample_dirs <- list.dirs(path = soupx_dir, full.names = TRUE, recursive = FALSE)
  print(soupx_sample_dirs)
  
  # Extract sample names and add pool ID 
  samples <- str_split_fixed(soupx_sample_dirs, "/", str_count(soupx_sample_dirs[1], "/")+1)[,str_count(soupx_sample_dirs[1], "/")+1]
  samples <- str_replace_all(samples, "-", ".")
  if (str_detect(pool, "pool")) {
    samples <- paste0(samples, "_", pool)
  }
  print(samples)
  
  # Create Seurat objects for all samples
  seurat_object_list <- list()
  for (j in 1:length(soupx_sample_dirs)) {
    seurat_object_list[[j]] <- load_seurat(soupx_sample_dirs[j], samples[j], TRUE)
  }
  
  # Merge objects
  s <- merge(seurat_object_list[[1]], unlist(seurat_object_list, use.names = FALSE)[2:length(seurat_object_list)],
             add.cell.ids = samples, project = "AN1792_scFRP")
  
  # Calculate percent mitochondrial expression 
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
  
  # Add pool and sample ID to meta data
  s$pool <- pool
  s$sample_id <- s$orig.ident
  s$sample_barcode <- row.names(s@meta.data)
  pool_list[[i]] <- s
  i <- i + 1
}

# Merge pool objects
s <- merge(pool_list[[1]], unlist(pool_list, use.names = FALSE)[2:length(pool_list)], project = "AN1792_scFRP")

# Save memory 
pool_list <- NULL
gc()

# Join and re-split layers based on unique sample ID 
s <- JoinLayers(s)
s[["RNA"]] <- split(s[["RNA"]], f = s$sample_id)

# Save pre-QC object
saveRDS(s, file = paste0(output_dir, "data/s_raw.rds"))

# Calculate general QC metrics per sample: strict min, flexible max for number of genes and UMI counts, strict MT max https://www.10xgenomics.com/analysis-guides/common-considerations-for-quality-control-filters-for-single-cell-rna-seq-data
# Only min UMI and nFeature thresholds are used to define low quality cells
thresholds <- data.frame()
for (sample in unique(s$sample_id)) {
  
  meta <- s@meta.data[s$sample_id == sample,]
  
  nfeat_median <- median(meta$nFeature_RNA)
  nfeat_min <- min(meta$nFeature_RNA)
  nfeat_min_thresh <- median(meta$nFeature_RNA) - 2*stats::mad(meta$nFeature_RNA) # 10X recommendation is 2-5x
  nfeat_max <- max(meta$nFeature_RNA)
  nfeat_max_thresh <- median(meta$nFeature_RNA) + 5*stats::mad(meta$nFeature_RNA)
  
  umi_median <- median(meta$nCount_RNA)
  umi_min <- min(meta$nCount_RNA)
  umi_min_thresh <- median(meta$nCount_RNA) - 3*stats::mad(meta$nCount_RNA) # 10X recommendation is 3-5x 
  umi_max <- max(meta$nCount_RNA)
  umi_max_thresh <- median(meta$nCount_RNA) + 5*stats::mad(meta$nCount_RNA)
  
  mt_median <- median(meta$percent.mt)
  mt_max <- max(meta$percent.mt)
  mt_max_thresh <- median(meta$percent.mt) + 3*stats::mad(meta$percent.mt) # 10X recommendation is 3-5x

  thresholds <- rbind(thresholds, data.frame(sample = sample, pool = unique(meta$pool), 
                                             umi_median = umi_median, umi_min = umi_min, 
                                             umi_min_thresh = umi_min_thresh,
                                             umi_max = umi_max, umi_max_thresh = umi_max_thresh, 
                                             nfeat_median = nfeat_median, nfeat_min = nfeat_min, nfeat_min_thresh = nfeat_min_thresh,
                                             nfeat_max = nfeat_max, nfeat_max_thresh = nfeat_max_thresh,
                                             mt_median = mt_median, mt_max = mt_max, mt_max_thresh = mt_max_thresh))
}

# Initialize table of QC stats 
qc_table <- data.frame(table(s$sample_id))
colnames(qc_table) <- c("sample", "pre_qc")
qc_table$pool <- s$pool[match(qc_table$sample, s$sample_id)]
qc_table <- qc_table[,c("sample", "pool", "pre_qc")]

# Filter cells using minimum UMI and nFeature thresholds
pre_umi <- table(s$sample_id)
cells_keep <- c()
for (sample in unique(s$sample_id)) {
  meta <- s@meta.data[s$sample_id == sample,]
  meta <- meta[meta$nFeature_RNA > thresholds$nfeat_min_thresh[thresholds$sample == sample] & meta$nCount_RNA > thresholds$umi_min_thresh[thresholds$sample == sample],]
  cells_keep <- c(cells_keep, rownames(meta))
}
s <- subset(s, cells = cells_keep)
umi_removed <- data.frame(pre_umi - table(s$sample_id))
qc_table$umi_nfeat_removed <- umi_removed$Freq

# Filter cells using uniform MT % threshold
pre_mt <- table(s$sample_id)
s <- subset(s, percent.mt < 20)
mt_removed <- data.frame(pre_mt - table(s$sample_id))
qc_table$mt_removed <- mt_removed$Freq

# Update QC summary table 
post_qc <- data.frame(table(s$sample_id))
qc_table$post_qc <- post_qc$Freq

# Save QC stats
write.csv(qc_table, paste0(output_dir, "data/qc_stats.csv"))

# Save object with low quality cells removed 
saveRDS(s, file = paste0(output_dir, "data/s_lowquality_removed.rds"))
