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
# Summary: Gene set enrichment analysis of microglia and macrophage DEGs with enrichR 
#
#-----------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("enrichR")
})

# Define output folder
output_folder <- "/path/to/all/regions/celltype/mast/output/"

# Define name of hallmark database 
dbs <- listEnrichrDbs()
hallmark_db <- "MSigDB_Hallmark_2020"  

# Load DE results
mac <- read.csv(paste0(output_folder, "results/Mac.csv"), row.names = 1)
mg <- read.csv(paste0(output_folder, "results/Mg.csv"), row.names = 1)

# Run enrichR on upregulated DEGs
mac_genes <- mac$gene[mac$DE == "Upregulated"]
mg_genes <- mg$gene[mg$DE == "Upregulated"]
mac_results <- enrichr(mac_genes, hallmark_db)
mg_results <- enrichr(mg_genes, hallmark_db)

# Filter for significant results 
mac_results <- mac_results[[1]]
mac_results <- mac_results[mac_results$Adjusted.P.value < 0.05,]
mg_results <- mg_results[[1]]
mg_results <- mg_results[mg_results$Adjusted.P.value < 0.05,]

# Save results
write.csv(mac_results, paste0(output_folder, "enrichR_upregulated/mac_all_regions.csv"))
write.csv(mg_results, paste0(output_folder, "enrichR_upregulated/mg_all_regions.csv"))

# Run enrichR on downregulated DEGs
mac_genes <- mac$gene[mac$DE == "Downregulated"]
mg_genes <- mg$gene[mg$DE == "Downregulated"]
mac_results <- enrichr(mac_genes, hallmark_db)
mg_results <- enrichr(mg_genes, hallmark_db)

# Filter for significant results 
mac_results <- mac_results[[1]]
mac_results <- mac_results[mac_results$Adjusted.P.value < 0.05,]
mg_results <- mg_results[[1]]
mg_results <- mg_results[mg_results$Adjusted.P.value < 0.05,]

# Save results
write.csv(mac_results, paste0(output_folder, "enrichR_downregulated/mac_all_regions.csv"))
write.csv(mg_results, paste0(output_folder, "enrichR_downregulated/mg_all_regions.csv"))

# Combine results for microglia and macrophages and identify top 10 overall upregulated and downregulated terms
mac1 <- read.csv(paste0(output_folder, "enrichR_downregulated/mac_all_regions.csv"), row.names = 1)
mg1 <- read.csv(paste0(output_folder, "enrichR_downregulated/mg_all_regions.csv"), row.names = 1)
mac1$type <- "Macrophages"
mg1$type <- "Microglia"
mac1$Combined.Score <- -mac1$Combined.Score
mg1$Combined.Score <- -mg1$Combined.Score
df1 <- rbind(mac1, mg1)

# Identify terms with top 10 overall scores 
df1 <- df1 %>% arrange(Combined.Score)
down_terms <- df1$Term[1:10]

mac2 <- read.csv(paste0(output_folder, "enrichR_upregulated/mac_all_regions.csv"), row.names = 1)
mg2 <- read.csv(paste0(output_folder, "enrichR_upregulated/mg_all_regions.csv"), row.names = 1)
mac2$type <- "Macrophages"
mg2$type <- "Microglia"
df2 <- rbind(mac2, mg2)

# Identify terms with top 10 overall scores 
df2 <- df2 %>% arrange(desc(Combined.Score))
up_terms <- df2$Term[1:10]

df <- rbind(df1, df2)
df <- df[df$Term %in% c(up_terms, down_terms),]

# Define cell type colors
colors <- c("#C77CFF", "#FDB863")
names(colors) <- c("Microglia", "Macrophages")

# Order results by combined score
df$plotorder <- rank(df$Combined.Score)
df <- df[order(df$plotorder),]

# Define variable for dot size
df$neg_logBH <- -log10(df$Adjusted.P.value)

# Generate dot plot
plt <- ggplot(df, aes(x = Combined.Score, y = reorder(Term, plotorder), fill = type, label = Term, size = neg_logBH)) +
  geom_point(stat="identity", stroke = NA, shape = 21) + 
  geom_text(size = 4.25, vjust = -100) +
  labs(fill = NULL) + ggtitle("enrichR MSigDB Analysis: LCMB vs. CAA (All Regions)") +
  theme_bw() + ylab("") + xlab("Combined Score") +
  scale_fill_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_size_continuous(range = c(3, 12)) + labs(size = "-log10(padj)") + 
  geom_vline(xintercept = 0, color = "black", linetype = "dotted")

pdf(paste0(output_folder, "enrichR_plots/combined_dotplot.pdf"), width = 8, height = 8)
print(plt)
dev.off()



