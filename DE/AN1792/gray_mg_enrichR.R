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
# Summary: Gene set enrichment analysis of DEGs in C2L microglia-enriched gray matter with enrichR 
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
output_folder <- "/path/to/gray/celltype/deseq2/output/folder/"

# Define name of hallmark database
dbs <- listEnrichrDbs()
hallmark_db <- "MSigDB_Hallmark_2020"  

# Filter operator 
`%notin%` <- Negate(`%in%`)

# Run enrichR for unique upregulated DEGs, unique downregulated DEGs, shared upregulated DEGs and shared downregulated DEGs

# Identify gene sets 

# Load results
lim <- read.csv(paste0(output_folder, "lim_vs_nAD/results/Microglia.csv"), row.names = 1)
ext <- read.csv(paste0(output_folder, "ext_vs_nAD/results/Microglia.csv"), row.names = 1)

# Remove results with NA adjusted p values 
lim <- lim[!is.na(lim$padj),]
ext <- ext[!is.na(ext$padj),]

# Subset for columns of interest
lim <- lim[,c("log2FoldChange", "padj", "gene", "DE")]
ext <- ext[,c("log2FoldChange", "padj", "gene", "DE")]

# Add data for genes that are not tested in both groups 
lim_add <- ext$gene[ext$gene %notin% lim$gene]
ext_add <- lim$gene[lim$gene %notin% ext$gene]

if (length(lim_add) > 0) {
  lim_add <- data.frame(log2FoldChange = rep(0, length(lim_add)), padj = rep(1, length(lim_add)),
                        gene = lim_add, DE = rep("Not DE", length(lim_add)), row.names = lim_add)
  lim <- rbind(lim, lim_add)
}

if (length(ext_add) > 0) {
  ext_add <- data.frame(log2FoldChange = rep(0, length(ext_add)), padj = rep(1, length(ext_add)),
                        gene = ext_add, DE = rep("Not DE", length(ext_add)), row.names = ext_add)
  ext <- rbind(ext, ext_add)
}

# Match gene order
lim <- lim[order(lim$gene),]
ext <- ext[order(ext$gene),]

# Calculate PFC 
lim$PFC <- -log10(lim$padj) * abs(lim$log2FoldChange)
ext$PFC <- -log10(ext$padj) * abs(ext$log2FoldChange)
print(sum(lim$PFC == Inf))
print(sum(ext$PFC == Inf))

# Combine results 
data <- data.frame(gene = lim$gene, lim_lfc = lim$log2FoldChange, ext_lfc = ext$log2FoldChange,
                   lim_de = lim$DE, ext_de = ext$DE, lim_pfc = lim$PFC, ext_pfc = ext$PFC,
                   row.names = lim$gene)

# Subset for genes that are DE in either or both groups 
data <- data[data$lim_de != "Not DE" | data$ext_de != "Not DE",]

# Create variable identifying significance in either or both groups (consider a gene "unique" if it changes direction)
data <- data %>% mutate(group = case_when(lim_de == "Upregulated" & ext_de == "Not DE"  ~ "Upregulated iAD-Lim", 
                                          lim_de == "Downregulated" & ext_de == "Not DE" ~ "Downregulated iAD-Lim", 
                                          ext_de == "Upregulated" & lim_de == "Not DE"  ~ "Upregulated iAD-Ext", 
                                          ext_de == "Downregulated" & lim_de == "Not DE" ~ "Downregulated iAD-Ext", 
                                          ext_de == "Upregulated" & lim_de == "Upregulated" ~ "Shared Upregulated",
                                          ext_de == "Downregulated" & lim_de == "Downregulated" ~ "Shared Downregulated",
                                          (ext_de != "Not DE" & lim_de != "Not DE") & (ext_de != lim_de) ~ "changes_direction")) 

# Update for DEGs that change direction (duplicate genes)
genes_to_add <- rownames(data)[data$group == "changes_direction"] 
if (length(genes_to_add) > 0) {
  lim <- data[genes_to_add,]
  rownames(lim) <- paste0("lim_", rownames(lim))
  ext <- data[genes_to_add,]
  rownames(ext) <- paste0("ext_", rownames(ext))
  lim$group <- ifelse(lim$lim_de == "Upregulated", "Upregulated iAD-Lim", "Downregulated iAD-Lim")
  ext$group <- ifelse(ext$ext_de == "Upregulated", "Upregulated iAD-Ext", "Downregulated iAD-Ext")
  data_add <- rbind(lim, ext)
  
  # Update data 
  data <- data[rownames(data) %notin% genes_to_add,]
  data <- rbind(data, data_add)
}

# Run enrichR for each group (manually update group: lim, ext, shared)
group_names <- c("Upregulated iAD-Lim", "Downregulated iAD-Lim", "Upregulated iAD-Ext", "Downregulated iAD-Ext",
                 "Shared Upregulated", "Shared Downregulated")
names(group_names) <- c("lim_up", "lim_down", "ext_up", "ext_down", "shared_up", "shared_down")
group <- "ext"

# Run enrichR for upregulated DEGs
up_degs <- data$gene[data$group == group_names[[paste0(group, "_up")]]]
up_terms <- enrichr(up_degs, hallmark_db)
up_terms <- up_terms[[1]]
up_terms <- up_terms[up_terms$Adjusted.P.value < 0.05,]
if (nrow(up_terms) > 0) {
  write.csv(up_terms, paste0(output_folder, "enrichR_shared_unique/microglia/", group, "_up.csv"))
}

# Run enrichR for downregulated DEGs
down_degs <- data$gene[data$group == group_names[[paste0(group, "_down")]]]
down_terms <- enrichr(down_degs, hallmark_db)
down_terms <- down_terms[[1]]
down_terms <- down_terms[down_terms$Adjusted.P.value < 0.05,]
if (nrow(down_terms) > 0) {
  write.csv(down_terms, paste0(output_folder, "enrichR_shared_unique/microglia/", group, "_down.csv"))
}


# Load results
group_names <- c("Unique iAD-Lim", "Unique iAD-Lim", "Unique iAD-Ext", "Unique iAD-Ext", "Shared", "Shared")
names(group_names) <- c("lim_up", "lim_down", "ext_up", "ext_down", "shared_up", "shared_down")
df1 <- data.frame()
df2 <- data.frame()
for (file in list.files(paste0(output_folder, "enrichR_shared_unique/microglia"))[!str_detect(list.files(paste0(output_folder, "enrichR_shared_unique/microglia")), ".pdf")]) {
  print(file)
  results <- read.csv(paste0(output_folder, "enrichR_shared_unique/microglia/", file), row.names = 1)
  
  # Negate combined score for analysis of downregulated genes
  if (str_detect(file, "down")) {
    results$group <- group_names[[str_replace(file, ".csv", "")]]
    results$Combined.Score <- -results$Combined.Score
    df2 <- rbind(df2, results)
  } else {
    results$group <- group_names[[str_replace(file, ".csv", "")]]
    df1 <- rbind(df1, results)
  }
}

# Combine results for upregulated and downregulated genes
df <- rbind(df1, df2)

# Define group colors
colors <- c("#C77CFF", "#FDB863", "#8DD3C7")
names(colors) <- c("Unique iAD-Lim", "Unique iAD-Ext", "Shared")

# Order results by combined score
df$plotorder <- rank(df$Combined.Score)
df <- df[order(df$plotorder),]

# Define variable for dot size
df$neg_logBH <- -log10(df$Adjusted.P.value)

# Generate dot plot
plt <- ggplot(df, aes(x = Combined.Score, y = reorder(Term, plotorder), fill = group, label = Term, size = neg_logBH)) +
  geom_point(stat="identity", stroke = NA, shape = 21) + 
  geom_text(size = 4.25, vjust = -100) +
  labs(fill = NULL) + ggtitle("enrichR MSigDB Analysis: Microglia-Enriched Gray Matter") +
  theme_bw() + ylab("") + xlab("Combined Score") +
  scale_fill_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) + scale_size_continuous(range = c(3, 12)) + labs(size = "-log10(padj)") + 
  geom_vline(xintercept = 0, color = "black", linetype = "dotted")

pdf(paste0(output_folder, "enrichR_shared_unique/microglia/combined_dotplot.pdf"), width = 8, height = 8)
print(plt)
dev.off()
