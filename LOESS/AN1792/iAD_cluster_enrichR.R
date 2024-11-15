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
# Summary: Gene set enrichment analysis of iAD LOESS clusters with enrichR 
#
#-----------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("enrichR")
  library("stringr")
})

# Define output folder
output_folder <- "/path/to/AN1792/loess/"

# Define name of hallmark database
dbs <- listEnrichrDbs()
hallmark_db <- "MSigDB_Hallmark_2020"  

# Load iAD cluster annotations
clusters <- read.csv(paste0(output_folder, "data/iAD_clusters_all_genes_broad.csv"), row.names = 1)
clusters$gene <- row.names(clusters)

# Run enrichR for each cluster (manually run cluster 1-12)
cluster <- 1
genes <- clusters$gene[clusters$nclust_12 == cluster]
terms <- enrichr(genes, hallmark_db)
terms <- terms[[1]]
terms <- terms[terms$Adjusted.P.value < 0.05,]
if (nrow(terms) > 0) {
  write.csv(terms, paste0(output_folder, "enrichR/iAD_all_genes/clust_", cluster, ".csv"))
}

# Load results
df <- data.frame()
for (file in list.files(paste0(output_folder, "enrichR/iAD_all_genes"))[!str_detect(list.files(paste0(output_folder, "enrichR/iAD_all_genes")), ".pdf")]) {
  print(file)
  results <- read.csv(paste0(output_folder, "enrichR/iAD_all_genes/", file), row.names = 1)
  results$clust <- str_replace(file, ".csv", "")
  df <- rbind(df, results)
}

# Define cluster colors
colors <- c("#8DD3C7", "#CC99FF", "#FDB462")
names(colors) <- paste0("clust_", c(1,3,4))

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

pdf(paste0(output_folder, "enrichR/iAD_all_genes/dotplot.pdf"), width = 8, height = 8)
print(plt)
dev.off()
