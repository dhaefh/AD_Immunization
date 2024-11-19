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
# Summary: Gene set enrichment analysis of lecanemab LOESS clusters with enrichR 
#
#-----------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("enrichR")
  library("stringr")
})

# Define output folder
output_folder <- "/path/to/lecanemab/loess/"

# Define name of hallmark database
dbs <- listEnrichrDbs()
hallmark_db <- "MSigDB_Hallmark_2020"  

# Load LCMB cluster annotations
clusters <- read.csv(paste0(output_folder, "data/LCMB_clusters_all_genes_broad.csv"), row.names = 1)
clusters$gene <- row.names(clusters)

# Run enrichR for 11 clusters 
for (cluster in 1:11) {
  genes <- clusters$gene[clusters$nclust_11 == cluster]
  print(paste0("Cluster ", cluster, ": ", length(genes), " genes"))
  terms <- enrichr(genes, hallmark_db)
  terms <- terms[[1]]
  terms <- terms[terms$Adjusted.P.value < 0.05,]
  print(nrow(terms))
  if (nrow(terms) > 0) {
    write.csv(terms, paste0(output_folder, "enrichR/clust_", cluster, ".csv"))
  }
  terms <- NULL
}


# Load results
df <- data.frame()
terms_keep <- c()
for (file in list.files(paste0(output_folder, "enrichR"))[!str_detect(list.files(paste0(output_folder, "enrichR")), ".pdf")]) {
  print(file)
  results <- read.csv(paste0(output_folder, "enrichR/", file), row.names = 1)
  results$clust <- str_replace(file, ".csv", "")
  
  # Keep top 13 terms if more than 13
  if (nrow(results) > 13) {
    results <- results %>% dplyr::arrange(desc(Combined.Score))
    terms <- results$Term[1:13]
  } else {
    terms <- results$Term
  }
  terms_keep <- c(terms_keep, terms)
  df <- rbind(df, results)
}

# Filter for top 13 terms per cluster
terms_keep <- unique(terms_keep)
df <- df[df$Term %in% terms_keep,]

# Define cluster colors
colors <- c("#8DD3C7", "#CC99FF", "#FDB462", "#3399CC")
names(colors) <- paste0("clust_", c(2,3,8,11))

# Order results by combined score
df$plotorder <- rank(df$Combined.Score)
df <- df[order(df$plotorder),]

# Define variable for dot size
df$neg_logBH <- -log10(df$Adjusted.P.value)

# Set factor levels for clusters
df$clust <- factor(df$clust, levels = names(colors))

# Generate dot plot
plt <- ggplot(df, aes(x = Combined.Score, y = reorder(Term, plotorder), fill = clust, label = Term, size = neg_logBH)) +
  geom_point(stat="identity", stroke = NA, shape = 21) + 
  geom_text(size = 4.25, vjust = -100) +
  labs(fill = NULL) + ggtitle("enrichR MSigDB Analysis: LOESS Clusters") +
  theme_bw() + ylab("") + xlab("Combined Score") +
  scale_fill_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) + scale_size_continuous(range = c(3, 12)) + labs(size = "-log10(padj)")

pdf(paste0(output_folder, "enrichR/dotplot.pdf"), width = 8, height = 8)
print(plt)
dev.off()
