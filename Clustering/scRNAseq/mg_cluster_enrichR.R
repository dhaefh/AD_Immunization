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
# Summary: Identify positive and negative markers for microglia clusters and perform gene set enrichment analysis with enrichR 
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("enrichR")
  library("stringr")
  library("Seurat")
})

# Define filter operator
`%notin%` <- Negate(`%in%`)

# Define output folder
output_dir <- "/path/to/microglia/enrichr/output/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load integrated scRNAseq immune Seurat object
s <- readRDS("/path/to/integrated/scRNAseq/immune/object.rds")

# Subset for lecanemab/nAD microglia clusters
s <- subset(s, immune_fine_updated %in% c("Mg-0", "Mg-1", "Mg-2", "Mg-3", "Mg-4") &
              sample_merged %notin% paste0("102.", c(1, 16, 17, 19, 22)))
gc()

# Recorrect SCT data
DefaultAssay(s) <- "SCT"
s <- PrepSCTFindMarkers(s)
Idents(s) <- "immune_fine_updated"

# Find positive and negative markers for each cluster
markers <- list()
for (cluster in unique(s$immune_fine_updated)){
  
  # Test genes expressed in 1% of both groups and 10% of either group (groups are current cluster vs. all other clusters)
  LFC <- FoldChange(s, ident.1 = cluster, slot = "data", assay = "SCT", base = 2)
  genes <- row.names(LFC)[(LFC$pct.1 >= 0.01 & LFC$pct.2 >= 0.01) & (LFC$pct.1 >= 0.1 | LFC$pct.2 >= 0.1)]
  
  # Exclude contamination genes
  if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes)) != 0) {
    genes <- genes[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes)] 
  }
  
  # Find markers
  results <- FindMarkers(s, ident.1 = cluster, features = genes, recorrect_umi = FALSE, min.pct = 0, only.pos = FALSE, 
                         fc.slot = "data", base = 2)
  results$BH <- p.adjust(results$p_val, method = "BH")
  results$cluster <- cluster
  results$gene <- rownames(results)
  markers[[cluster]] <- results
}

# Compile and save results
all_markers <- data.table::rbindlist(markers)|>as.data.frame()
write.csv(all_markers, paste0(output_dir, "markers.csv"), row.names = FALSE, quote = FALSE)

# Identify name of hallmark database
dbs <- listEnrichrDbs()
hallmark_db <- "MSigDB_Hallmark_2020"  

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Run enrichR for positive markers
for (clust in unique(all_markers$cluster)) {
  print(clust)
  clust_markers <- all_markers[all_markers$cluster == clust,]
  genes <- clust_markers$gene[clust_markers$BH < p_thresh & clust_markers$avg_log2FC > fc_thresh]
  print(length(genes))
  terms <- enrichr(genes, hallmark_db)
  terms <- terms[[1]]
  terms <- terms[terms$Adjusted.P.value < 0.05,]
  if (nrow(terms) > 0) {
    write.csv(terms, paste0(output_dir, clust, "_up.csv"))
  }
  terms <- NULL
  gc()
}

# Run enrichR for negative markers
for (clust in unique(all_markers$cluster)) {
  print(clust)
  clust_markers <- all_markers[all_markers$cluster == clust,]
  genes <- clust_markers$gene[clust_markers$BH < p_thresh & clust_markers$avg_log2FC < -fc_thresh]
  print(length(genes))
  terms <- enrichr(genes, hallmark_db)
  terms <- terms[[1]]
  terms <- terms[terms$Adjusted.P.value < 0.05,]
  if (nrow(terms) > 0) {
    write.csv(terms, paste0(output_dir, clust, "_down.csv"))
  }
  terms <- NULL
  gc()
}

# Load and combine results, filtering for top 5 upregulated and downregulated terms per cluster
df <- data.frame()
terms_keep <- c()
for (file in list.files(output_dir)[!str_detect(list.files(output_dir), ".pdf")]) {
  if (file == "markers.csv") {
    next
  }
  print(file)
  results <- read.csv(paste0(output_dir, file), row.names = 1)
  results$cluster <- str_split_fixed(file, "_", 2)[,1]
  
  if (str_detect(file, "_up")) {
    results <- results %>% dplyr::arrange(desc(Combined.Score))
    if (nrow(results) > 5) {
      terms <- results$Term[1:5]
    } else {
      terms <- results$Term
    }
  } else {
    # Negate combined score for analysis of negative markers
    results$Combined.Score <- -results$Combined.Score 
    results <- results %>% dplyr::arrange(Combined.Score)
    if (nrow(results) > 5) {
      terms <- results$Term[1:5]
    } else {
      terms <- results$Term
    }
  }
  terms_keep <- c(terms_keep, terms)
  df <- rbind(df, results)
}
df <- df[df$Term %in% unique(terms_keep),]

# Define cluster colors
colors <- c("#8DD3C7", "#E78AC3", "#CC99FF", "#FDB462", "#3399CC")
names(colors) <- paste0("Mg-", 0:4)
df$cluster <- factor(df$cluster, levels = names(colors))

# Order results by combined score
df$plotorder <- rank(df$Combined.Score)
df <- df[order(df$plotorder),]

# Define variable for dot size
df$neg_logBH <- -log10(df$Adjusted.P.value)

# Generate dot plot
plt <- ggplot(df, aes(x = Combined.Score, y = reorder(Term, plotorder), fill = cluster, label = Term, size = neg_logBH)) +
  geom_point(stat="identity", stroke = NA, shape = 21) + 
  geom_text(size = 4.25, vjust = -100) +
  labs(fill = NULL) + ggtitle("enrichR MSigDB Analysis: Microglia Clusters") +
  theme_bw() + ylab("") + xlab("Combined Score") +
  scale_fill_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) + scale_size_continuous(range = c(3, 12)) + labs(size = "-log10(padj)") + 
  geom_vline(xintercept = 0, color = "black", linetype = "dotted")

pdf(paste0(output_dir, "combined_dotplot.pdf"), width = 12, height = 8)
print(plt)
dev.off()