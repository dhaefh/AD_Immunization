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
  library("extrafont")
  library("showtext")
})

# Define output folder
output_folder <- "/path/to/manual/layer/deseq2/output/folder/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Load integrated AN1792 Seurat object
s <- readRDS("/path/to/integrated/AN1792/object.rds")

# Combine counts across samples
DefaultAssay(s) <- "Spatial"
s <- JoinLayers(s)

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Define "clusters" (manual layers) and comparisons to run 
clusters <- unique(s$manual_annotation)
comparisons <- c("iAD_vs_nAD", "nAD_vs_NNC", "ext_vs_nAD", "lim_vs_nAD")

# Run DESeq2 for each manual layer for each comparison
for (comparison in comparisons) {
  
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
  
  for (cluster in clusters) {
    
    # Subset for manual layer and set default assay 
    clust_s <- subset(s, manual_annotation == cluster)
    DefaultAssay(clust_s) <- "Spatial"
    
    # Extract raw counts
    counts <- GetAssayData(clust_s, assay = "Spatial", layer = "counts")
    
    # If comparing iAD-Ext or iAD-Lim, set idents to condition_residual_amyloid
    if ((ident.1 %in% c("ext", "lim")) | (ident.2 %in% c("ext", "lim"))) {
      Idents(clust_s) <- "condition_residual_amyloid"
    } else {
      Idents(clust_s) <- "condition"
    }
    
    # Calculate percent expression for comparison groups based on raw counts
    LFC <- FoldChange(object = clust_s, ident.1 = ident.1, ident.2 = ident.2, fc.name = "log2FoldChange", base = 2, 
                      assay = "Spatial", slot = "counts") %>% data.frame()
    
    # Filter for genes expressed in 1% of either group
    genes_keep <- row.names(LFC)[LFC$pct.1 >= 0.01 | LFC$pct.2 >= 0.01]
    
    # Exclude contamination genes
    if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)) != 0) {
      genes_keep <- genes_keep[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)] 
    }
    
    # Filter counts for genes to test
    counts <- counts[genes_keep,]
    
    # Extract meta data for manual layer
    meta <- s@meta.data[s$manual_annotation == cluster,]
    
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
    meta <- meta[,c('age', 'sex', 'sample_id', 'avg_feat', 'gDNA_percent', 'condition', 'condition_residual_amyloid')]
    meta <- unique(meta)
    meta$age_centered <- scale(meta$age)
    meta$feat_centered <- scale(meta$avg_feat)
    meta$gDNA_centered <- scale(meta$gDNA_percent)
    meta$sex <- as.factor(meta$sex)
    row.names(meta) <- meta$sample_id
    
    # Ensure pseudobulk data and meta data have the same sample order
    bulk <- bulk[, row.names(meta)]
    print(sum(colnames(bulk) != row.names(meta)))
    print(sum(is.na(bulk)))
    
    # Update condition variable if comparing iAD-Ext or iAD-Lim
    if ((ident.1 %in% c("ext", "lim")) | (ident.2 %in% c("ext", "lim"))) {
      meta$condition <- factor(meta$condition_residual_amyloid)
    } else {
      meta$condition <- factor(meta$condition)
    }
    
    # Set condition reference level
    meta$condition <- relevel(meta$condition, ref = ident.2)
    
    # Run DESeq2 
    skip_to_next <- FALSE
    tryCatch(
      {
        dds <- DESeqDataSetFromMatrix(countData = bulk, colData = meta, design= ~ sex + age_centered + feat_centered + gDNA_centered + condition)
        
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
      print(paste0("Error in ", comparison, ", ", cluster))
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
}




