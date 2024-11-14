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
# Summary: Run SCTransform and perform CCA integration of immune cells in Seurat
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
})

# Load integrated scRNAseq Seurat object
s <- readRDS("/path/to/integrated/scRNAseq/object.rds")

# Subset for initially defined immune subtypes
s <- subset(s, broad_anno_long %in% c("microglia/macrophage", "T lymphocyte", "monocyte"))
DefaultAssay(s) <- "RNA"
s <- DietSeurat(s, assays = "RNA")
gc()

# Create variable for integration groups (AN1792, CAA, LCMB)
s$pool_merged <- as.character(s$pool)
s$pool_merged[grep("^AN1792", s$pool)] <- "AN1792"
s$pool_merged[grep("^LCMB", s$pool)] <- "LCMB"
unique(s$pool_merged)

# Split layers by merged sample ID 
s <- JoinLayers(s)
s[["RNA"]] <- split(s[["RNA"]], f = s$pool_merged)

# Re-run SCTransform, PCA and IntegrateLayers on immune subset
s <- SCTransform(s, assay = "RNA", vars.to.regress = "percent.mt")
s <- RunPCA(s, verbose = F)
s <- IntegrateLayers(object = s, method = CCAIntegration, normalization.method = 'SCT', orig.reduction = 'pca', verbose = F, new.reduction = "cca")

# Remove old clustering from meta data
s@meta.data[,grep("harmony", colnames(s@meta.data))] <- NULL
s$seurat_clusters <- NULL

# Generate UMAPs with 20 and 30 PCs (30 PCs are ultimately used for final UMAP + clustering)
ElbowPlot(s, ndims = 50)
s <- RunUMAP(s, dims = 1:20, reduction = "cca", reduction.name = "ccaumap20pcs")
s <- RunUMAP(s, dims = 1:30, reduction = "cca", reduction.name = "ccaumap30pcs")

# Save object 
saveRDS(s, "/path/to/integrated/scRNAseq/immune/object.rds")




