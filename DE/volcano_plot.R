# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                         AN1792 Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Written by: Anne Forsyth 
# Summary: Generate volcano plot using DESeq2 or MAST output
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ggplot2")
  library("ggrepel")
  library("UpSetR")
  library("randomcoloR")
})

# Load results
results <- read.csv("/path/to/DE/results.csv", row.names = 1)

# Specify column names corresponding to LFC, adjusted p-value and gene name
lfc <- "<name of LFC column>"
padj <- "<name of padj column>"
gene <- "<name of gene name column>"

# Define DE thresholds
fc_thresh <- log2(1.5) 
p_thresh <- 0.05

# Define plot function 
volcano <- function(results, title = NULL, label_genes = NULL) {
  
  # Define plot colors
  colors <- c("blue", "red", "black")
  names(colors) <- c("Downregulated", "Upregulated", "Not DE")
  
  # Rename variables 
  results$lfc <- results[[lfc]]
  results$padj <- results[[padj]]
  results$gene <- results[[gene]]
  
  # Calculate PFC 
  results$neg_logBH <- -log10(results$padj)
  results$neg_logBH[results$neg_logBH == Inf] <- max(results$neg_logBH[results$neg_logBH != Inf]) + 10
  results$PFC <- results$neg_logBH*abs(results$lfc)
  
  # Define DEGs
  results$DE <- "Not DE"
  results$DE[results$padj < p_thresh & results$lfc > fc_thresh] <- "Upregulated"
  results$DE[results$padj < p_thresh & results$lfc < -fc_thresh] <- "Downregulated"
  
  # Gene labeling 
  results$label <- NA
  
  # If genes not specified, label top 30 by PFC, otherwise label specified genes
  if (is.null(label_genes)) {
    results <- results %>% dplyr::arrange(desc(PFC))
    degs <- results$gene[results$DE != "Not DE"][1:30]
    results$label[results$gene %in% degs] <- results$gene[results$gene %in% degs]
  } else {
    results$label[results$gene %in% label_genes] <- results$gene[results$gene %in% label_genes]
  }
  
  # Define basic title if not specified 
  if (is.null(title)) {
    title <- paste0(sum(results$DE != "Not DE"), " DEGs")
  }
  
  # Generate plot 
  plt <- ggplot(data.frame(results), aes(x = lfc, y = neg_logBH, fill = DE, label = label)) + 
    geom_point(aes(size = PFC), alpha = 0.5, shape = 21, stroke = NA) + 
    scale_size_continuous(range = c(3, 12)) + 
    geom_label_repel(color = "black", size = 6, force = 3, min.segment.length = 0, force_pull = 0, box.padding = 1, max.overlaps = Inf,
                     fill = "white", alpha = 0.75, label.size = NA) +
    geom_vline(xintercept=c(-fc_thresh, fc_thresh), linetype = 3) + geom_hline(yintercept= -log10(p_thresh), linetype = 3) +
    scale_fill_manual(values = colors) + 
    theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black'), aspect.ratio = 1) + 
    xlim(min(results$lfc[!is.na(results$padj)]) - 0.2, max(results$lfc[!is.na(results$padj)]) + 0.2) +
    ggtitle(paste0(title)) + theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x = "Log2 Fold Change", y = "-log10(Adjusted P-Value)") +
    theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
  
  return(plt)
}
