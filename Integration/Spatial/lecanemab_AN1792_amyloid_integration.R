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
# Summary: Run SCTransform and perform Harmony integration in Seurat for lecanemab and AN1792 amyloid-rich spots in gray matter
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages ({
  library('Seurat')
  library('dplyr')
})

# Filter operator
`%notin%` <- Negate(`%in%`)

# Load all-cohort integrated object 
s <- readRDS("/path/to/integrated/all/cohorts/object.rds")

# Load cohort 5/7/8 data
cohort578 <- readRDS("/path/to/cohort578/integrated/object.rds")
cohort578 <- cohort578@meta.data
gc()

# Load cohort 1 data
cohort1 <- readRDS("/path/to/cohort1/integrated/object.rds")
cohort1 <- cohort1@meta.data
gc()

# Define condition variable for cohort 5/7/8
cohort578$condition <- NA
cohort578$condition[grep("^A14|^A11|^NMA22.A", cohort578$sample_id)] <- "CAA"
cohort578$condition[grep("^NMA22.B", cohort578$sample_id)] <- "LCMB"
print(unique(cohort578[,c("sample_id", "condition")]))

# Define gray matter amyloid enrichment for both groups
cohort578$gray_amyloid <- "not_rich"
cohort578$gray_amyloid[cohort578$cortical_amyloid > 183 & cohort578$vascular_amyloid == 0] <- "rich" # All spots with cortical amyloid density > 183 and vascular amyloid density = 0
cohort578$gray_amyloid[-grep("^gray", cohort578$manual_layer)] <- "not_rich" # Exclude anything not in gray matter
cohort1$gray_amyloid <- "not_rich"
cohort1$gray_amyloid[cohort1$amyloid_filter == "include" & cohort1$amyloid_fluo > 183 & cohort1$vessel_neighbor == "not_vessel"] <- "rich" # Amyloid "include", amyloid density > 183, exclude vascular amyloid + neighbors
cohort1$gray_amyloid[-grep("^gray", cohort1$manual_annotation)] <- "not_rich" # Exclude anything not in gray matter

# Combine meta data
add_meta <- data.frame(condition = c(cohort1$condition, cohort578$condition), condition_clearance = c(cohort1$condition_clearance, cohort578$condition),
                       gray_amyloid = c(cohort1$gray_amyloid, cohort578$gray_amyloid), row.names = c(rownames(cohort1), rownames(cohort578)))
add_meta <- add_meta[rownames(s@meta.data),]
s@meta.data <- cbind(s@meta.data, add_meta)

# Subset for gray matter amyloid-rich spots, remove 2 donors with only 1 spot (error in SCTransform)
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

# Save object
saveRDS(s, "/path/to/integrated/all/cohorts/amyloid/object.rds")


