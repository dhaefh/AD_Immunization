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
# Summary: Remove doublets and create post-QC Seurat object for lecanemab
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
output_dir <- "/path/to/lecanemab/QC/output/"

# Load sample objects with doublet annotations
s_list <- list()
i <- 1
objects <- list.files(paste0(output_dir, "data/DoubletFinder/"))
for (object in objects) {
  s_list[[i]] <- readRDS(paste0(output_dir, "data/DoubletFinder/", object))
  i <- i + 1
}

# Manually combine raw counts and meta data from all samples (more efficient than merge() for number of samples)
s <- s_list[[1]]
names(s@assays$RNA@layers)[grep("^counts", names(s@assays$RNA@layers))] <- "counts"
counts <- GetAssayData(s, assay = "RNA", layer = "counts")
meta <- s@meta.data
colnames(meta)[grep("^pANN_0.25", colnames(meta))] <- "pANN_0.25"
meta[,c("RNA_snn_res.0.8", "seurat_clusters")] <- NULL
for (i in 2:length(s_list)) {
  s <- s_list[[i]]
  names(s@assays$RNA@layers)[grep("^counts", names(s@assays$RNA@layers))] <- "counts"
  cur_counts <- GetAssayData(s, assay = "RNA", layer = "counts")
  print(sum(rownames(cur_counts) != rownames(counts)))
  counts <- cbind(counts, cur_counts)
  
  cur_meta <- s@meta.data
  colnames(cur_meta)[grep("^pANN_0.25", colnames(cur_meta))] <- "pANN_0.25"
  cur_meta[,c("RNA_snn_res.0.8", "seurat_clusters")] <- NULL
  print(sum(colnames(cur_meta) != colnames(meta)))
  meta <- rbind(meta, cur_meta)
}

# Create Seurat object with merged counts
s <- CreateSeuratObject(counts = counts, assay = "RNA", project = "AN1792_scFRP")

# Save memory
s_list <- NULL
gc()

# Update meta data
sum(meta$nFeature_RNA != s$nFeature_RNA)
sum(meta$nCount_RNA != s$nCount_RNA)
sum(rownames(meta) != rownames(s@meta.data))
s$orig.ident <- meta$orig.ident
s@meta.data <- cbind(s@meta.data, meta[,4:ncol(meta)])
sum(meta != s@meta.data)
Idents(s) <- "orig.ident"

# Split layers by sample ID
s[["RNA"]] <- split(s[["RNA"]], f = s$sample_id)

# Extract meta data with full annotations before removing doublets (not adjusted for homotypic proportion)
meta <- s@meta.data
s <- subset(s, DF.unadj == "Singlet")

# Save post-QC object
saveRDS(s, file = paste0(output_dir, "data/s_post_qc.rds"))

# Save doublet summary data
write.csv(table(meta$sample_id, meta$DF.unadj), paste0(output_dir, "data/doublet_table.csv"))
write.csv(meta[,c("sample_id", "pool", "sample_barcode", "DF.unadj", "DF.adj")],
          paste0(output_dir, "data/doublet_annotations.csv"))


