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

# Load libraries
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
})

# Define input and output folders
input_folder <- "/path/to/preprocessing/output/folder/"
output_folder <- "/path/to/qc/output/folder/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_folder, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_folder, "plots"), showWarnings = FALSE, recursive = TRUE)

# Load binned Seurat objects for each sample 
s_list <- list()
i <- 1
for (sample in list.files(input_folder)) {
  
  # Load sample object
  s <- readRDS(paste0(input_folder, sample, "/binned_seurat.rds"))
  
  # Add nuclei areas and sample ID to meta data
  areas <- read.csv(paste0(input_folder, sample, "/areas.csv"))
  rownames(areas) <- areas$id
  areas <- areas[rownames(s@meta.data),]
  sum(is.na(areas)) %>% print()
  sum(rownames(areas) != rownames(s@meta.data)) %>% print()
  s$area <- areas$area
  s$sample_id <- sample
  
  # Add to list
  s_list[[i]] <- s
  i <- i + 1
}

# Merge objects
s <- merge(s_list[[1]], unlist(s_list, use.names = FALSE)[2:length(s_list)], project = "AN1792_HIPP")

# Re-split layers based on sample ID
s <- JoinLayers(s)
s[['RNA']] <- split(s[['RNA']], f = s$sample_id)

# QC plots
qc_table <- data.frame(table(s$sample_id))
colnames(qc_table) <- c("sample_id", "pre_qc")
rownames(qc_table) <- qc_table$sample_id

# Keep nuclei with UMI > 10 
pre_umi <- table(s$sample_id)
s <- subset(s, nCount_RNA > 10)
umi_removed <- data.frame(pre_umi - table(s$sample_id))
qc_table$umi_removed <- umi_removed$Freq

# Calculate percent mitochondrial expression 
s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")

# Keep nuclei with MT % < 20 
pre_mt <- table(s$sample_id)
s <- subset(s, percent.mt < 20)
mt_removed <- data.frame(pre_mt - table(s$sample_id))
qc_table$mt_removed <- mt_removed$Freq

# Save QC stats
post_qc <- data.frame(table(s$sample_id))
qc_table$post_qc <- post_qc$Freq
qc_table$sample_id <- NULL
write.csv(qc_table, paste0(output_folder, "data/qc_stats.csv"))

# Save post-QC object
saveRDS(s, paste0(output_folder, "data/s_post_qc.rds"))