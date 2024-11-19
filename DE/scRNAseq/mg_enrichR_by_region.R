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
# Summary: Gene set enrichment analysis of microglia DEGs with enrichR 
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("enrichR")
  library("stringr")
})

# Define output folder
output_folder <- "/path/to/by/region/celltype/mast/output/"

# Define name of hallmark database 
dbs <- listEnrichrDbs()
hallmark_db <- "MSigDB_Hallmark_2020"  

# Filter operator 
`%notin%` <- Negate(`%in%`)

# Run enrichR for upregulated DEGs in each brain region 
for (comp in list.files(output_folder)[!str_detect(list.files(output_folder), ".csv") & list.files(output_folder) %notin% c("hallmark", "enrichR_upregulated", "enrichR_downregulated")]) {
  results <- read.csv(paste0(output_folder, comp, "/results/Mg.csv"), row.names = 1)
  genes <- results$gene[results$DE == "Upregulated"]
  results <- enrichr(genes, hallmark_db)
  results <- results[[1]]
  results <- results[results$Adjusted.P.value < 0.05,]
  write.csv(results, paste0(output_folder, "enrichR_upregulated/", comp, ".csv"))
}

# Run enrichR for downregulated DEGs in each brain region
for (comp in list.files(output_folder)[!str_detect(list.files(output_folder), ".csv") & list.files(output_folder) %notin% c("hallmark", "enrichR_upregulated", "enrichR_downregulated")]) {
  results <- read.csv(paste0(output_folder, comp, "/results/Mg.csv"), row.names = 1)
  genes <- results$gene[results$DE == "Downregulated"]
  results <- enrichr(genes, hallmark_db)
  results <- results[[1]]
  results <- results[results$Adjusted.P.value < 0.05,]
  write.csv(results, paste0(output_folder, "enrichR_downregulated/", comp, ".csv"))
}

# Load downregulated results
df1 <- data.frame()
for (file in list.files(paste0(output_folder, "enrichR_downregulated"))[!str_detect(list.files(paste0(output_folder, "enrichR_downregulated")), ".pdf")]) {
  results <- read.csv(paste0(output_folder, "enrichR_downregulated/", file), row.names = 1)
  if (nrow(results) > 0) {
    results$comp <- str_replace(file, ".csv", "")
    df1 <- rbind(df1, results)
  }
}

# Negate combined score
df1$Combined.Score <- -df1$Combined.Score

# Load upregulated results
df2 <- data.frame()
for (file in list.files(paste0(output_folder, "enrichR_upregulated"))[!str_detect(list.files(paste0(output_folder, "enrichR_downregulated")), ".pdf")]) {
  results <- read.csv(paste0(output_folder, "enrichR_upregulated/", file), row.names = 1)
  if (nrow(results) > 0) {
    results$comp <- str_replace(file, ".csv", "")
    df2 <- rbind(df2, results)
  }
}

# Filter for top 10 overall upregulated terms, keep all downregulated terms (less than 10 total) 
df2 <- df2 %>% arrange(desc(Combined.Score))
up_terms <- df2$Term[1:10]
down_terms <- df1$Term

df <- rbind(df1, df2)
df <- df[df$Term %in% c(up_terms, down_terms),]

# Define region colors
colors <- c("#C77CFF", "#FDB863", "#8DD3C7", "#E78AC3")
names(colors) <- c("B1_vs_A1", "B3_vs_A3", "B4_vs_A4", "B9_vs_A9")

# Order results by enrichment score
df$plotorder <- rank(df$Combined.Score)
df <- df[order(df$plotorder),]

# Define variable for dot size
df$neg_logBH <- -log10(df$Adjusted.P.value)

# Generate dot plot
plt <- ggplot(df, aes(x = Combined.Score, y = reorder(Term, plotorder), fill = comp, label = Term, size = neg_logBH)) +
  geom_point(stat="identity", stroke = NA, shape = 21) + 
  geom_text(size = 4.25, vjust = -100) +
  labs(fill = NULL) + ggtitle("enrichR MSigDB Analysis: Microglia") +
  theme_bw() + ylab("") + xlab("Combined Score") +
  scale_fill_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) + scale_size_continuous(range = c(3, 12)) + labs(size = "-log10(padj)") + 
  geom_vline(xintercept = 0, color = "black", linetype = "dotted")

pdf(paste0(output_folder, "enrichr_combined_dotplot.pdf"), width = 8, height = 8)
print(plt)
dev.off()
