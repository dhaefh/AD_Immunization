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
# Summary: Partial Spearman correlation between iAD microglia-associated gene expression and clinical hallmarks 
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages ({
  library(Seurat)
  library(dplyr)
  library(plyr)
  library(data.table)
  library(tidyverse)
  library(Seurat)
  library(rlist)
  library(ggrepel)
  library(ggplot2)
  library(ggpubr)
  library(ppcor)
})

# Define output folder 
output_dir <- "/path/to/clinical/correlation/output/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define filter operator
`%notin%` <- Negate(`%in%`)

# Load integrated AN1792 Seurat object
s <- readRDS("/path/to/integrated/AN1792/object.rds")

# Filter for iAD microglia-enriched ST spots in gray matter 
s <- subset(s, manual_annotation %notin% c("white", "meninges") & Microglia_enriched == 1 & condition == "iAD")
gc()

# Combine counts across samples
DefaultAssay(s) <- "Spatial"
s <- JoinLayers(s)

# Calculate average log-normalized TREM2 and APOE expression per sample 
data <- GetAssayData(s, assay = "Spatial", layer = "data")
data <- data[c("TREM2", "APOE"),]
data <- as.matrix(data)
sum(colnames(data) != row.names(s@meta.data))
colnames(data) <- s$sample_id
aggregated <- t(rowsum(t(data), group = colnames(data)))
for (col in colnames(aggregated)) {
  aggregated[,col] <- aggregated[,col] / sum(colnames(data) == col)
}
aggregated <- t(aggregated)

# Extract meta data and calculate average nFeatures per sample
meta <- s@meta.data
nfeat_summary <- meta %>% dplyr::group_by(sample_id) %>% dplyr::summarize(sample_nfeat = mean(nFeature_Spatial))
meta$avg_feat <- NA
for (sample in unique(meta$sample_id)) {
  meta$avg_feat[meta$sample_id == sample] <- nfeat_summary$sample_nfeat[nfeat_summary$sample_id == sample]
}

# Subset metadata to have one row per sample 
meta <- meta[,c('sample_id', 'age', 'sex', 'avg_feat', 'gDNA_percent')]
meta <- unique(meta)
rownames(meta) <- meta$sample_id

# Add clinical hallmark meta data 
hallmarks <- read.csv("/path/to/cohort1_clinical_data.csv", row.names = 1)
hallmarks <- hallmarks[,3:ncol(hallmarks)]
colnames(hallmarks) <- c("mean_titre", "survival_time", "dementia_duration", "plaque_score", "pct_amy_gray", "pct_tau_gray")
rownames(hallmarks) <- str_replace_all(rownames(hallmarks), "-", ".")
rownames(hallmarks) <- str_replace_all(rownames(hallmarks), "_", ".")
hallmarks <- hallmarks[rownames(meta),]
meta <- cbind(meta, hallmarks)

# Standardize continuous covariates and make sex numeric
meta$age_centered <- scale(meta$age)
meta$feat_centered <- scale(meta$avg_feat)
meta$gDNA_centered <- scale(meta$gDNA_percent)
meta$sex_numeric <- ifelse(meta$sex == "f", 1, 0)

# Add gene expression to meta data 
aggregated <- aggregated[rownames(meta),]
meta <- cbind(meta, aggregated)

# Run Spearman correlation for each clinical hallmark 
for (hallmark in c("mean_titre", "survival_time", "dementia_duration", "plaque_score", "pct_amy_gray", "pct_tau_gray")) {
  
  # Initialize data frame for results
  results <- data.frame()
  for (gene in c("TREM2", "APOE")) {
    
    # Partial Spearman correlation
    ctest <- ppcor::pcor.test(x = meta[,hallmark], 
                              y = meta[,gene], 
                              z = meta[,c("sex_numeric", "age_centered", "feat_centered", "gDNA_centered")],
                              method = "spearman")
    
    # Update results
    results <- rbind(results, data.frame(gene = gene, spearmans_rho = ctest$estimate, p = ctest$p.value, row.names = gene))
  }
  write.csv(results, paste0(output_dir, hallmark, ".csv"))
}

# Define significant correlations
rho_thresh <- 0.1
p_thresh <- 0.05
for (file in list.files(output_dir)) {
  print(file)
  results <- read.csv(paste0(output_dir, file), row.names = 1)
  results$padj <- p.adjust(results$p, method = "BH")
  results$direction <- "Not Significant"
  results$direction[results$padj < p_thresh & results$spearmans_rho > rho_thresh] <- "Positive"
  results$direction[results$padj < p_thresh & results$spearmans_rho < -rho_thresh] <- "Negative"
  write.csv(results, paste0(output_dir, file))
}



