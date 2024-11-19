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
# Summary: Summarize amyloid and microglia associated genes for AN1792 and Lecanemab based on ranking within DE analyses
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
})

# Define output folder
output_dir <- "/path/to/gene/ranking/output/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load DE results

# Amyloid-rich gray matter
an1792_gm_plaque <- read.csv("/path/to/cohort1/amyloid/mast/output/folder/iAD_vs_nAD/results/amyloid_rich.csv", row.names = 1)
lcmb_gm_plaque <- read.csv("/path/to/cohort578/amyloid/all/regions/mast/output/folder/results/amyloid_rich.csv", row.names = 1)

# Plaque cluster 6 
an1792_plaque6 <- read.csv("/path/to/cluster6/mast/output/iAD_vs_nAD/results/plaque_6.csv", row.names = 1)
lcmb_plaque6 <- read.csv("/path/to/cluster6/mast/output/LCMB_vs_CAA/results/plaque_6.csv", row.names = 1)

# Microglia
an1792_mg <- read.csv("/path/to/cohort1/general/celltype/deseq2/output/folder/iAD_vs_nAD/results/Microglia.csv", row.names = 1)
lcmb_mg <- read.csv("/path/to/scRNAseq/all/regions/celltype/mast/output/results/Mg.csv", row.names = 1)

# C2L enriched microglia clusters in plaque niche
an1792_mg2 <- read.csv("/path/to/mg/cluster/mast/output/lim_vs_nAD/results/Mg_2.csv", row.names = 1)
an1792_mg4 <- read.csv("/path/to/mg/cluster/mast/output/lim_vs_nAD/results/Mg_4.csv", row.names = 1)
lcmb_mg2 <- read.csv("/path/to/mg/cluster/mast/output/LCMB_vs_CAA/results/Mg_2.csv", row.names = 1)
lcmb_mg4 <- read.csv("/path/to/mg/cluster/mast/output/LCMB_vs_CAA/results/Mg_4.csv", row.names = 1)

# AN1792 microglia 
summary <- data.frame()
mg_results <- list(an1792_mg, an1792_mg2, an1792_mg4)
names(mg_results) <- c("AN1792 Mg", "AN1792 Mg-2", "AN1792 Mg-4")
for (index in names(mg_results)) {
  results <- mg_results[[index]]
  
  # Calculate signed PFC and rank genes (ensuring PFC values are finite)
  if ("log2FoldChange" %in% colnames(results)) {
    results <- results[!is.na(results$padj),]
    results$neg_logBH <- -log10(results$padj)
    results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
    print(sum(results$neg_logBH == Inf))
    results$PFC <- results$neg_logBH*results$log2FoldChange
  } else {
    results <- results[!is.na(results$BH),]
    results$neg_logBH <- -log10(results$BH)
    results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
    print(sum(results$neg_logBH == Inf))
    results$PFC <- results$neg_logBH*results$avg_log2FC
  }
  
  # Rank genes by PFC, calculate percentile, filter for positive fold change
  results <- results %>% dplyr::arrange(PFC)
  results$rank <- 1:nrow(results)
  results$percentile <- results$rank/nrow(results)
  results <- results[results$PFC > 0,]
  results$de <- ifelse(results$DE == "Not DE", 0, 1)
  
  # Update summary table
  summary <- rbind(summary, data.frame(gene = results$gene, PFC = results$PFC, DE = results$de,
                                       rank = results$rank, percentile = results$percentile, analysis = index))
}

# Save full results 
write.csv(summary, paste0(output_dir, "mg_an1792_summary.csv"), row.names = FALSE)

# Summarize number of occurences and average percentile of each gene 
gene_summary <- summary %>% dplyr::group_by(gene) %>% dplyr::summarize(n = n(), n_de = sum(DE), avg_rank = mean(percentile))
gene_summary <- gene_summary %>% dplyr::arrange(desc(avg_rank))
write.csv(gene_summary, paste0(output_dir, "mg_an1792_ranking.csv"), row.names = FALSE)

# LCMB microglia
summary <- data.frame()
mg_results <- list(lcmb_mg, lcmb_mg2, lcmb_mg4)
names(mg_results) <- c("LCMB Mg", "LCMB Mg-2", "LCMB Mg-4")
for (index in names(mg_results)) {
  results <- mg_results[[index]]
  
  # Calculate signed PFC and rank genes (ensuring PFC values are finite)
  if ("log2FoldChange" %in% colnames(results)) {
    results <- results[!is.na(results$padj),]
    results$neg_logBH <- -log10(results$padj)
    results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
    print(sum(results$neg_logBH == Inf))
    results$PFC <- results$neg_logBH*results$log2FoldChange
  } else {
    results <- results[!is.na(results$BH),]
    results$neg_logBH <- -log10(results$BH)
    results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
    print(sum(results$neg_logBH == Inf))
    results$PFC <- results$neg_logBH*results$avg_log2FC
  }
  
  # Rank genes by PFC, calculate percentile, filter for positive fold change
  results <- results %>% dplyr::arrange(PFC)
  results$rank <- 1:nrow(results)
  results$percentile <- results$rank/nrow(results)
  results <- results[results$PFC > 0,]
  results$de <- ifelse(results$DE == "Not DE", 0, 1)
  
  # Update summary table
  summary <- rbind(summary, data.frame(gene = results$gene, PFC = results$PFC, DE = results$de,
                                       rank = results$rank, percentile = results$percentile, analysis = index))
}

# Save full results 
write.csv(summary, paste0(output_dir, "mg_lcmb_summary.csv"), row.names = FALSE)

# Summarize number of occurences and average percentile of each gene 
gene_summary <- summary %>% dplyr::group_by(gene) %>% dplyr::summarize(n = n(), n_de = sum(DE), avg_rank = mean(percentile))
gene_summary <- gene_summary %>% dplyr::arrange(desc(avg_rank))
write.csv(gene_summary, paste0(output_dir, "mg_lcmb_ranking.csv"), row.names = FALSE)

# AN1792 amyloid
summary <- data.frame()
amyloid_results <- list(an1792_gm_plaque, an1792_plaque6)
names(amyloid_results) <- c("AN1792 Plaque", "AN1792 Plaque-6")
for (index in names(amyloid_results)) {
  results <- amyloid_results[[index]]
  
  # Calculate signed PFC and rank genes (ensuring PFC values are finite)
  if ("log2FoldChange" %in% colnames(results)) {
    results <- results[!is.na(results$padj),]
    results$neg_logBH <- -log10(results$padj)
    results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
    print(sum(results$neg_logBH == Inf))
    results$PFC <- results$neg_logBH*results$log2FoldChange
  } else {
    results <- results[!is.na(results$BH),]
    results$neg_logBH <- -log10(results$BH)
    results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
    print(sum(results$neg_logBH == Inf))
    results$PFC <- results$neg_logBH*results$avg_log2FC
  }
  
  # Rank genes by PFC, calculate percentile, filter for positive fold change
  results <- results %>% dplyr::arrange(PFC)
  results$rank <- 1:nrow(results)
  results$percentile <- results$rank/nrow(results)
  results <- results[results$PFC > 0,]
  results$de <- ifelse(results$DE == "Not DE", 0, 1)
  
  # Update summary table
  summary <- rbind(summary, data.frame(gene = results$gene, PFC = results$PFC, DE = results$de,
                                       rank = results$rank, percentile = results$percentile, analysis = index))
}

# Save full results 
write.csv(summary, paste0(output_dir, "amyloid_an1792_summary.csv"), row.names = FALSE)

# Summarize number of occurences and average percentile of each gene 
gene_summary <- summary %>% dplyr::group_by(gene) %>% dplyr::summarize(n = n(), n_de = sum(DE), avg_rank = mean(percentile))
gene_summary <- gene_summary %>% dplyr::arrange(desc(avg_rank))
write.csv(gene_summary, paste0(output_dir, "amyloid_an1792_ranking.csv"), row.names = FALSE)

# LCMB amyloid
summary <- data.frame()
amyloid_results <- list(lcmb_gm_plaque, lcmb_plaque6)
names(amyloid_results) <- c("LCMB Plaque", "LCMB Plaque-6")
for (index in names(amyloid_results)) {
  results <- amyloid_results[[index]]
  
  # Calculate signed PFC and rank genes (ensuring PFC values are finite)
  if ("log2FoldChange" %in% colnames(results)) {
    results <- results[!is.na(results$padj),]
    results$neg_logBH <- -log10(results$padj)
    results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
    print(sum(results$neg_logBH == Inf))
    results$PFC <- results$neg_logBH*results$log2FoldChange
  } else {
    results <- results[!is.na(results$BH),]
    results$neg_logBH <- -log10(results$BH)
    results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
    print(sum(results$neg_logBH == Inf))
    results$PFC <- results$neg_logBH*results$avg_log2FC
  }
  
  # Rank genes by PFC, calculate percentile, filter for positive fold change
  results <- results %>% dplyr::arrange(PFC)
  results$rank <- 1:nrow(results)
  results$percentile <- results$rank/nrow(results)
  results <- results[results$PFC > 0,]
  results$de <- ifelse(results$DE == "Not DE", 0, 1)
  
  # Update summary table
  summary <- rbind(summary, data.frame(gene = results$gene, PFC = results$PFC, DE = results$de,
                                       rank = results$rank, percentile = results$percentile, analysis = index))
}

# Save full results 
write.csv(summary, paste0(output_dir, "amyloid_lcmb_summary.csv"), row.names = FALSE)

# Summarize number of occurences and average percentile of each gene 
gene_summary <- summary %>% dplyr::group_by(gene) %>% dplyr::summarize(n = n(), n_de = sum(DE), avg_rank = mean(percentile))
gene_summary <- gene_summary %>% dplyr::arrange(desc(avg_rank))
write.csv(gene_summary, paste0(output_dir, "amyloid_lcmb_ranking.csv"), row.names = FALSE)

# AN1792 microglia + amyloid
summary <- data.frame()
combined_results <- list(an1792_mg, an1792_mg2, an1792_mg4, an1792_gm_plaque, an1792_plaque6)
names(combined_results) <- c("AN1792 Mg", "AN1792 Mg-2", "AN1792 Mg-4", "AN1792 Plaque", "AN1792 Plaque-6")
for (index in names(combined_results)) {
  results <- combined_results[[index]]
  
  # Calculate signed PFC and rank genes (ensuring PFC values are finite)
  if ("log2FoldChange" %in% colnames(results)) {
    results <- results[!is.na(results$padj),]
    results$neg_logBH <- -log10(results$padj)
    results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
    print(sum(results$neg_logBH == Inf))
    results$PFC <- results$neg_logBH*results$log2FoldChange
  } else {
    results <- results[!is.na(results$BH),]
    results$neg_logBH <- -log10(results$BH)
    results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
    print(sum(results$neg_logBH == Inf))
    results$PFC <- results$neg_logBH*results$avg_log2FC
  }
  
  # Rank genes by PFC, calculate percentile, filter for positive fold change
  results <- results %>% dplyr::arrange(PFC)
  results$rank <- 1:nrow(results)
  results$percentile <- results$rank/nrow(results)
  results <- results[results$PFC > 0,]
  results$de <- ifelse(results$DE == "Not DE", 0, 1)
  
  # Update summary table
  summary <- rbind(summary, data.frame(gene = results$gene, PFC = results$PFC, DE = results$de,
                                       rank = results$rank, percentile = results$percentile, analysis = index))
}

# Save full results 
write.csv(summary, paste0(output_dir, "combined_an1792_summary.csv"), row.names = FALSE)

# Summarize number of occurences and average percentile of each gene 
gene_summary <- summary %>% dplyr::group_by(gene) %>% dplyr::summarize(n = n(), n_de = sum(DE), avg_rank = mean(percentile))
gene_summary <- gene_summary %>% dplyr::arrange(desc(avg_rank))
write.csv(gene_summary, paste0(output_dir, "combined_an1792_ranking.csv"), row.names = FALSE)

# LCMB microglia + amyloid
summary <- data.frame()
combined_results <- list(lcmb_mg, lcmb_mg2, lcmb_mg4, lcmb_gm_plaque, lcmb_plaque6)
names(combined_results) <- c("LCMB Mg", "LCMB Mg-2", "LCMB Mg-4", "LCMB Plaque", "LCMB Plaque-6")
for (index in names(combined_results)) {
  results <- combined_results[[index]]
  
  # Calculate signed PFC and rank genes (ensuring PFC values are finite)
  if ("log2FoldChange" %in% colnames(results)) {
    results <- results[!is.na(results$padj),]
    results$neg_logBH <- -log10(results$padj)
    results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
    print(sum(results$neg_logBH == Inf))
    results$PFC <- results$neg_logBH*results$log2FoldChange
  } else {
    results <- results[!is.na(results$BH),]
    results$neg_logBH <- -log10(results$BH)
    results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
    print(sum(results$neg_logBH == Inf))
    results$PFC <- results$neg_logBH*results$avg_log2FC
  }
  
  # Rank genes by PFC, calculate percentile, filter for positive fold change
  results <- results %>% dplyr::arrange(PFC)
  results$rank <- 1:nrow(results)
  results$percentile <- results$rank/nrow(results)
  results <- results[results$PFC > 0,]
  results$de <- ifelse(results$DE == "Not DE", 0, 1)
  
  # Update summary table
  summary <- rbind(summary, data.frame(gene = results$gene, PFC = results$PFC, DE = results$de,
                                       rank = results$rank, percentile = results$percentile, analysis = index))
}

# Save full results 
write.csv(summary, paste0(output_dir, "combined_lcmb_summary.csv"), row.names = FALSE)

# Summarize number of occurences and average percentile of each gene 
gene_summary <- summary %>% dplyr::group_by(gene) %>% dplyr::summarize(n = n(), n_de = sum(DE), avg_rank = mean(percentile))
gene_summary <- gene_summary %>% dplyr::arrange(desc(avg_rank))
write.csv(gene_summary, paste0(output_dir, "combined_lcmb_ranking.csv"), row.names = FALSE)
