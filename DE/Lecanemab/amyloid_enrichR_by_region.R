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
# Summary: Gene set enrichment analysis of DEGs in amyloid-rich gray matter with enrichR 
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
output_folder <- "/path/to/amyloid/by/region/mast/output/folder/"

# Define name of hallmark database
dbs <- listEnrichrDbs()
hallmark_db <- "MSigDB_Hallmark_2020"  

# Filter operator 
`%notin%` <- Negate(`%in%`)

# Run enrichR for each comparison

# Manually change comparison (B1_vs_A1, B3_vs_A3, B4_vs_A4, B9_vs_A9)
comp <- "B9_vs_A9"

# Run enrichR for upregulated DEGs
results <- read.csv(paste0(output_folder, comp, "/results/amyloid_rich.csv"), row.names = 1)
up_degs <- results$gene[results$DE == "Upregulated"]
up_terms <- enrichr(up_degs, hallmark_db)
up_terms <- up_terms[[1]]
up_terms <- up_terms[up_terms$Adjusted.P.value < 0.05,]
if (nrow(up_terms) > 0) {
  write.csv(up_terms, paste0(output_folder, "enrichR/", comp, "_up.csv"))
}

# Run enrichR for downregulated DEGs
down_degs <- results$gene[results$DE == "Downregulated"]
down_terms <- enrichr(down_degs, hallmark_db)
down_terms <- down_terms[[1]]
down_terms <- down_terms[down_terms$Adjusted.P.value < 0.05,]
if (nrow(down_terms) > 0) {
  write.csv(down_terms, paste0(output_folder, "enrichR/", comp, "_down.csv"))
}


# Load results
df1 <- data.frame()
df2 <- data.frame()
for (file in list.files(paste0(output_folder, "enrichR"))[!str_detect(list.files(paste0(output_folder, "enrichR")), ".pdf")]) {
  print(file)
  results <- read.csv(paste0(output_folder, "enrichR/", file), row.names = 1)
  
  # Negate combined score for analysis of downregulated genes
  if (str_detect(file, "down")) {
    results$comp <- str_replace(file, "_down.csv", "")
    results$Combined.Score <- -results$Combined.Score
    df2 <- rbind(df2, results)
  } else {
    results$comp <- str_replace(file, "_up.csv", "")
    df1 <- rbind(df1, results)
  }
}

# Filter for top 10 unique upregulated terms for B3 vs. A3, keep all terms for other regions - note that there are no significant downregulated terms
`%notin%` <- Negate(`%in%`)
table(df1$comp)
temp <- df1[df1$comp == "B3_vs_A3",]
temp <- temp %>% dplyr::arrange(desc(Combined.Score))
terms <- temp$Term[temp$Term %notin% df1$Term[df1$comp != "B3_vs_A3"]][1:10] # top 10 unique terms for B3 vs. A3
terms <- c(terms, temp$Term[temp$Term %in% df1$Term[df1$comp != "B3_vs_A3"]]) # keep terms shared with other comparisons
df1 <- df1[df1$comp != "B3_vs_A3" | df1$Term %in% terms,] # only data removed is unique to B3 vs. A3
df <- rbind(df1, df2)

# Dot colors
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
  labs(fill = NULL) + ggtitle("enrichR MSigDB Analysis: Amyloid-Rich Gray Matter") +
  theme_bw() + ylab("") + xlab("Combined Score") +
  scale_fill_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) + scale_size_continuous(range = c(3, 12)) + labs(size = "-log10(padj)") 

pdf(paste0(output_folder, "enrichR/combined_dotplot.pdf"), width = 8, height = 8)
print(plt)
dev.off()
