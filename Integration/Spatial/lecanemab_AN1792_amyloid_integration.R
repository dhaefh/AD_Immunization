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
# Summary: Run SCTransform and perform Harmony integration in Seurat for lecanemab and AN1792 Aß-rich spots in gray matter
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages ({
  library('Seurat')
  library('dplyr')
})

# Filter operator
`%notin%` <- Negate(`%in%`)

# Load all-cohort integrated object 
s <- readRDS("/path/to/integrated/all/cohorts/object.rds")

# Load lecanemab data
cohort578 <- readRDS("/path/to/lecanemab/integrated/object.rds")
cohort578 <- cohort578@meta.data
gc()

# Load AN1792 data
cohort1 <- readRDS("/path/to/AN1792/integrated/object.rds")
cohort1 <- cohort1@meta.data
gc()

# Define gray matter Aß enrichment for both groups
cohort578$gray_amyloid <- "not_rich"
cohort578$gray_amyloid[cohort578$amyloid_neighbor_final == "amyloid"] <- "rich" 
cohort578$gray_amyloid[-grep("^gray", cohort578$manual_layer)] <- "not_rich"
cohort1$gray_amyloid <- "not_rich"
cohort1$gray_amyloid[cohort1$amyloid_neighbor_final == "amyloid"] <- "rich"
cohort1$gray_amyloid[-grep("^gray", cohort1$manual_annotation)] <- "not_rich"

# Combine meta data
add_meta <- data.frame(condition = c(cohort1$condition, cohort578$condition), condition_clearance = c(cohort1$condition_clearance, cohort578$condition),
                       gray_amyloid = c(cohort1$gray_amyloid, cohort578$gray_amyloid), row.names = c(rownames(cohort1), rownames(cohort578)))
add_meta <- add_meta[rownames(s@meta.data),]
s@meta.data <- cbind(s@meta.data, add_meta)

# Subset for gray matter Aß-rich spots, remove 2 donors with only 1 spot (error in SCTransform)
s <- subset(s, gray_amyloid == "rich" & condition %in% c("CAA", "LCMB", "iAD", "nAD") & sample_id %notin% c("AN1792.102.8", "AN1792.102.16"))
gc()
DefaultAssay(s) <- "Spatial"
s <- DietSeurat(s, assays = "Spatial")
s <- JoinLayers(s)
s[["Spatial"]] <- split(s[["Spatial"]], f = s$sample_id)

# Run SCTransform
s <- SCTransform(s, assay = 'Spatial')

# Run PCA
s <- RunPCA(s, verbose = F)

# Integrate data 
s <- IntegrateLayers(object = s, method = HarmonyIntegration, normalization.method = 'SCT', orig.reduction = 'pca', verbose = F, new.reduction = "harmony")

# Generate UMAP 
ElbowPlot(s, ndims = 50)
s <- RunUMAP(s, dims = 1:30, reduction = "harmony", reduction.name = "harmonyumap")

# Save object
saveRDS(s, "/path/to/integrated/all/cohorts/amyloid/object.rds")


