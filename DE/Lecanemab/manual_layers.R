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
# Summary: Differential expression with DESeq2 for manual layers
#
#-----------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("DESeq2")
  library("ggrepel")
  library("UpSetR")
  library("MAST")
})

# Define output folder
output_folder <- "/path/to/manual/layer/deseq2/output/folder/"

# Load integrated cohort 5/7/8 Seurat object
s <- readRDS("/path/to/integrated/cohort578/object.rds")

# Define filter operator 
`%notin%` <- Negate(`%in%`)

# Define gDNA %
s$gDNA_percent <- NA
s$gDNA_percent[s$sample_id == "NMA22.A1"] <- 1.4
s$gDNA_percent[s$sample_id == "NMA22.A3"] <- 1.3
s$gDNA_percent[s$sample_id == "NMA22.A4"] <- 2.5
s$gDNA_percent[s$sample_id == "NMA22.A9"] <- 1.8
s$gDNA_percent[s$sample_id == "NMA22.B1"] <- 1.0
s$gDNA_percent[s$sample_id == "NMA22.B3"] <- 0.1
s$gDNA_percent[s$sample_id == "NMA22.B4"] <- 0.2
s$gDNA_percent[s$sample_id == "NMA22.B9"] <- 0.5
s$gDNA_percent[s$sample_id == "A14.193.1"] <- 15.3
s$gDNA_percent[s$sample_id == "A14.193.3"] <- 16.3
s$gDNA_percent[s$sample_id == "A14.193.4"] <- 10.1
s$gDNA_percent[s$sample_id == "A14.193.9"] <- 1.8
s$gDNA_percent[s$sample_id == "A11.170.1"] <- 4.4
s$gDNA_percent[s$sample_id == "A11.170.3"] <- 4.9
s$gDNA_percent[s$sample_id == "A11.170.4"] <- 0.6
s$gDNA_percent[s$sample_id == "A11.170.9"] <- 0.8
print(sum(is.na(s$gDNA_percent)))

# Define region 
s$region <- NA
s$region[s$sample_id %in% paste0(c("A14.193.", "A11.170.", "NMA22.A", "NMA22.B"), 1)] <- "FCX"
s$region[s$sample_id %in% paste0(c("A14.193.", "A11.170.", "NMA22.A", "NMA22.B"), 3)] <- "TCX"
s$region[s$sample_id %in% paste0(c("A14.193.", "A11.170.", "NMA22.A", "NMA22.B"), 4)] <- "PCX"
s$region[s$sample_id %in% paste0(c("A14.193.", "A11.170.", "NMA22.A", "NMA22.B"), 9)] <- "HIPP"

# Subset for non-hippocampus
s <- subset(s, region %in% c("FCX", "TCX", "PCX"))
gc()

# Define condition variable
s$condition <- NA
s$condition[grep("\\.B", s$sample_id)] <- "LCMB"
s$condition[grep("\\.A|^A", s$sample_id)] <- "CAA"

# Combine counts across samples
DefaultAssay(s) <- "Spatial"
s <- JoinLayers(s)

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Define name of DE comparison
comparison <- "LCMB_vs_CAA"

# Extract names of comparison groups
ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]

# Define comparison name for DESeq2
comp_name <- paste0("condition_", comparison)

# Create subfolders
dir.create(paste0(output_folder, comparison), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_folder, comparison, "/results"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_folder, comparison, "/dispersion_plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_folder, comparison, "/plots"), showWarnings = FALSE, recursive = TRUE)

# Run DESeq2 for each manual layer
for (cluster in unique(s$manual_layer)) {
  
  # Subset for manual layer and set default assay
  clust_s <- subset(s, manual_layer == cluster)
  DefaultAssay(clust_s) <- "Spatial"
  
  # Extract raw counts
  counts <- GetAssayData(clust_s, assay = "Spatial", layer = "counts")
  
  # Set idents
  Idents(clust_s) <- "condition"
  
  # Calculate percent expression for comparison groups based on raw counts
  LFC <- FoldChange(object = clust_s, ident.1 = ident.1, ident.2 = ident.2, fc.name = "log2FoldChange", base = 2, 
                    assay = "Spatial", layer = "counts") %>% data.frame()
  
  # Filter for genes expressed in 1% of either group
  genes_keep <- row.names(LFC)[LFC$pct.1 >= 0.01 | LFC$pct.2 >= 0.01]
  
  # Exclude contamination genes
  if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)) != 0) {
    genes_keep <- genes_keep[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)] 
  }
  
  # Filter counts for genes to test
  counts <- counts[genes_keep,]
  
  # Extract meta data for manual layer
  meta <- s@meta.data[s$manual_layer == cluster,]
  
  # Calculate average nFeatures per sample within the current manual layer
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
  meta <- meta[,c('sample_id', 'region', 'avg_feat', 'gDNA_percent', 'condition')]
  meta <- unique(meta)
  meta$feat_centered <- scale(meta$avg_feat)
  meta$gDNA_centered <- scale(meta$gDNA_percent)
  meta$region <- factor(meta$region)
  row.names(meta) <- meta$sample_id
  
  # Ensure pseudobulk data and meta data have the same sample order
  bulk <- bulk[, row.names(meta)]
  print(sum(colnames(bulk) != row.names(meta)))
  print(sum(is.na(bulk)))
  
  # Set condition reference level
  meta$condition <- factor(meta$condition)
  meta$condition <- relevel(meta$condition, ref = ident.2)
  
  # Run DESeq2 
  skip_to_next <- FALSE
  tryCatch(
    {
      dds <- DESeqDataSetFromMatrix(countData = bulk, colData = meta, design= ~ region + feat_centered + gDNA_centered + condition)
      
      # Generate plots for local vs. parametric fit, and calculate median absolute residuals 
      parametric <- DESeq(dds)
      pdf(paste0(output_folder, comparison, "/dispersion_plots/", cluster, "_parametric.pdf"), height = 10, width = 10)
      plotDispEsts(parametric)
      dev.off()
      parametric_residual <- median(abs(mcols(parametric)$dispGeneEst - mcols(parametric)$dispFit))
      
      local <- DESeq(dds, fitType = "local")
      pdf(paste0(output_folder, comparison, "/dispersion_plots/", cluster, "_local.pdf"), height = 10, width = 10)
      plotDispEsts(local)
      dev.off()
      local_residual <- median(abs(mcols(local)$dispGeneEst - mcols(local)$dispFit))
      
      write.csv(data.frame(local = local_residual, parametric = parametric_residual), paste0(output_folder, comparison, "/dispersion_plots/", cluster, "_median_absolute_residual.csv"), row.names = FALSE)
      
      # Run DESeq with local fit
      dds <- DESeq(dds, fitType = "local")
    },
    error = function(e) {
      skip_to_next <<- TRUE
    }
  )
  if (skip_to_next) {
    print(paste0("Error in ", cluster))
    next
  }
  
  tryCatch(
    {
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
      write.csv(results, paste0(output_folder, comparison, "/results/", cluster, ".csv"))
    },
    error = function(e) {
      print(paste0("No ", comparison, " results for ", cluster))
    }
  )
}


