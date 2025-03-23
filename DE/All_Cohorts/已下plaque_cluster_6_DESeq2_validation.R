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
# Summary: Differential expression with DESeq2 for Aß cluster 6
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("DESeq2")
  library("ggrepel")
  library("UpSetR")
  library("randomcoloR")
  library("ComplexUpset")
})

# Define output folder
output_folder <- "/path/to/plaque6/deseq2/output/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Load all-cohorts integrated Seurat object and subset for cluster 6
s <- readRDS("/path/to/integrated/all/cohorts/object.rds")
temp <- readRDS("/path/to/integrated/all/cohorts/amyloid/object.rds")
spots_keep <- rownames(temp@meta.data[temp$plaque_cluster == "Aß-6",])
s <- subset(s, cells = spots_keep)
temp <- NULL
gc()

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Combine counts across samples
DefaultAssay(s) <- "Spatial"
s <- JoinLayers(s)

# Define comparisons to run 
comparisons <- c("iAD_vs_nAD", "LCMB_vs_CAA")

# Run DESeq2 
for (comparison in comparisons) {
  
  if (comparison == "iAD_vs_nAD") {
    cur_s <- subset(s, condition %in% c("iAD", "nAD"))
    print(unique(cur_s$condition))
  } else {
    cur_s <- subset(s, condition %in% c("LCMB", "CAA"))
    print(unique(cur_s$condition))
  }
  
  # Extract names of comparison groups
  ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
  ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
  
  # Define comparison name for DESeq2
  comp_name <- paste0("condition_", comparison)
  
  # Create subfolders
  dir.create(paste0(output_folder, comparison), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(output_folder, comparison, "/results"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(output_folder, comparison, "/plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(output_folder, comparison, "/dispersion_plots"), showWarnings = FALSE, recursive = TRUE)
  
  # Extract raw counts
  counts <- GetAssayData(cur_s, assay = "Spatial", layer = "counts")
  
  # Set idents 
  Idents(cur_s) <- "condition"
  
  # Calculate percent expression for comparison groups based on raw counts
  LFC <- FoldChange(object = cur_s, ident.1 = ident.1, ident.2 = ident.2, fc.name = "avg_log2FC", base = 2, 
                    assay = "Spatial", slot = "counts") %>% data.frame()
  
  # Filter for genes expressed in 1% of either group
  genes_keep <- row.names(LFC)[LFC$pct.1 >= 0.01 | LFC$pct.2 >= 0.01]
  
  # Exclude contamination genes
  if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)) != 0) {
    genes_keep <- genes_keep[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)] 
  }
  
  # Filter counts for genes to test
  counts <- counts[genes_keep,]
  
  # Extract meta data 
  meta <- cur_s@meta.data
  
  # Calculate average nFeatures per sample
  nfeat_summary <- meta %>% group_by(sample_id) %>% summarize(sample_nfeat = mean(nFeature_Spatial))
  meta$avg_feat <- NA
  for (sample in unique(meta$sample_id)) {
    meta$avg_feat[meta$sample_id == sample] <- nfeat_summary$sample_nfeat[nfeat_summary$sample_id == sample]
  }
  print(sum(is.na(meta$avg_feat)))
  
  # Confirm rows of meta data match columns of counts
  print(sum(colnames(counts) != row.names(meta)))
  
  # Sum raw counts per sample to generate pseudobulk data
  colnames(counts) <- meta$sample_id
  bulk <- t(rowsum(t(counts), group = colnames(counts)))
  
  # Subset metadata to have one row per sample and standardize continuous covariates
  meta <- meta[,c('sample_id', 'age', 'sex', 'avg_feat', 'gDNA_percent', 'region', 'condition')]
  meta <- unique(meta)
  meta$age_centered <- scale(meta$age)
  meta$feat_centered <- scale(meta$avg_feat)
  meta$sex <- factor(meta$sex)
  meta$gDNA_centered <- scale(meta$gDNA_percent)
  meta$region <- factor(meta$region)
  row.names(meta) <- meta$sample_id
  
  # Ensure pseudobulk data and meta data have the same sample order
  bulk <- bulk[, row.names(meta)]
  print(sum(colnames(bulk) != row.names(meta)))
  print(sum(is.na(bulk)))
  
  # Set reference level for condition
  meta$condition <- factor(meta$condition)
  meta$condition <- relevel(meta$condition, ref = ident.2)
  
  # Run DESeq2 with cohort-specific covariates
  if (comparison == "iAD_vs_nAD") {
    print(paste0("Including age and sex, excluding region for ", comparison))
    dds <- DESeqDataSetFromMatrix(countData = bulk, colData = meta, design= ~ sex + age_centered + feat_centered + gDNA_centered + condition)
  } else {
    print(paste0("Excluding age and sex, including region for ", comparison))
    dds <- DESeqDataSetFromMatrix(countData = bulk, colData = meta, design= ~ feat_centered + region + gDNA_centered + condition)
  }
  
  # Generate plots for local vs. parametric fit, and calculate median absolute residuals 
  parametric <- DESeq(dds)
  pdf(paste0(output_folder, comparison, "/dispersion_plots/parametric.pdf"), height = 10, width = 10)
  plotDispEsts(parametric)
  dev.off()
  parametric_residual <- median(abs(mcols(parametric)$dispGeneEst - mcols(parametric)$dispFit))
  
  local <- DESeq(dds, fitType = "local")
  pdf(paste0(output_folder, comparison, "/dispersion_plots/local.pdf"), height = 10, width = 10)
  plotDispEsts(local)
  dev.off()
  local_residual <- median(abs(mcols(local)$dispGeneEst - mcols(local)$dispFit))
  
  write.csv(data.frame(local = local_residual, parametric = parametric_residual), paste0(output_folder, comparison, "/dispersion_plots/median_absolute_residual.csv"), row.names = FALSE)
  
  # Run DESeq with local fit
  dds <- DESeq(dds, fitType = "local")
  
  # Generate results without independent filtering, using default outlier filtering
  results <- results(dds, name = paste0(comp_name), independentFiltering = FALSE)
  
  # Apply apeglm LFC shrinkage 
  results <- lfcShrink(dds, coef=paste0(comp_name), type="apeglm", res = results)
  
  # Define DEGs and save results
  results$gene <- row.names(results)
  results$DE <- "Not DE"
  results$DE[results$log2FoldChange > fc_thresh & results$padj < p_thresh] <- "Upregulated"
  results$DE[results$log2FoldChange < -fc_thresh & results$padj < p_thresh] <- "Downregulated"
  results$DE_gene <- NA
  results$DE_gene[results$DE != "Not DE"] <- results$gene[results$DE != "Not DE"]
  write.csv(results, paste0(output_folder, comparison, "/results/plaque_6.csv"))
}

# Define list of genes to validate and re-adjust p-values
validation_genes <- c("SPP1", "TREM2", "RAB13", "TSPO", "FAM107A", "S100A4", "APOE", "A2M", "CHI3L1", "TMSB4X", "FCGBP", "CTSB") 

# iAD vs. nAD
results <- read.csv(paste0(output_folder, "iAD_vs_nAD/results/plaque_6.csv"), row.names = 1)
results <- results[validation_genes,]
results$padj_filtered <- p.adjust(results$pvalue, method = "BH")
write.csv(results, paste0(output_folder, "iAD_vs_nAD/results/plaque_6_readjusted.csv"))

# LCMB vs. nAD
results <- read.csv(paste0(output_folder, "LCMB_vs_CAA/results/plaque_6.csv"), row.names = 1)
results <- results[validation_genes,]
results$padj_filtered <- p.adjust(results$pvalue, method = "BH")
write.csv(results, paste0(output_folder, "LCMB_vs_CAA/results/plaque_6_readjusted.csv"))
