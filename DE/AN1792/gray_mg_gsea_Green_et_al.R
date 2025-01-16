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
# Summary: Microglia state enrichment analysis of C2L microglia-enriched gray matter with fgsea, using human microglia states from Green et al. 
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library(readxl)
  library(Seurat)
  library(SeuratObject)
  library(tidyverse)
  library(cowplot)
  library(fgsea)
})

# Define output folder 
output_folder <- "/path/to/Green_et_al/gsea/output/folder/"

# Load microglia state markers
markers <- read_excel("/path/to/Green_et_al/states.xlsx")

# Filter operator 
`%notin%` <- Negate(`%in%`)

# Filter for microglia states 
markers <- markers[grep("^Mic", markers$state),]

# Group markers by microglia state
markers <- group_by(markers, state)

# Load DE results
de_folder <- "/path/to/gray/celltype/deseq2/output/folder/"
ext_nAD <- read.csv(paste0(de_folder, "ext_vs_nAD/results/Microglia.csv"))
iAD_nAD <- read.csv(paste0(de_folder, "iAD_vs_nAD/results/Microglia.csv"))
nAD_NNC <- read.csv(paste0(de_folder, "nAD_vs_NNC/results/Microglia.csv"))
lim_nAD <- read.csv(paste0(de_folder, "lim_vs_nAD/results/Microglia.csv"))
iAD_nAD <- iAD_nAD[!is.na(iAD_nAD$padj),]
ext_nAD <- ext_nAD[!is.na(ext_nAD$padj),]
nAD_NNC <- nAD_NNC[!is.na(nAD_NNC$padj),]
lim_nAD <- lim_nAD[!is.na(lim_nAD$padj),]

# Calculate signed PFC and ensure values are finite
iAD_nAD$PFC <- -log10(iAD_nAD$padj)*iAD_nAD$log2FoldChange
nAD_NNC$PFC <- -log10(nAD_NNC$padj)*nAD_NNC$log2FoldChange
ext_nAD$PFC <- -log10(ext_nAD$padj)*ext_nAD$log2FoldChange
lim_nAD$PFC <- -log10(lim_nAD$padj)*lim_nAD$log2FoldChange
print(sum(nAD_NNC$PFC == Inf))
print(sum(ext_nAD$PFC == Inf))
print(sum(lim_nAD$PFC == Inf))
print(sum(iAD_nAD$PFC == Inf))

# Get list of state markers
marklist <- lapply(split(select(markers, gene), markers$state), deframe)
names(marklist)

# Compile list of comparison results
complist <- list(nAD_NNC, ext_nAD, lim_nAD, iAD_nAD)
names(complist) <- c("nAD_NNC", "ext_nAD", "lim_nAD", "iAD_nAD")

# Format title names
titles <- c("nAD vs. NNC", "iAD-Ext vs. nAD", "iAD-Lim vs. nAD", "iAD vs. nAD")
names(titles) <- names(complist)

# Run GSEA for each comparison
options(mc.cores = 1)
gsealist <- list()
for (comp in names(complist)) {
  ranked_degs <- deframe(select(complist[[comp]], gene, PFC)) # Get gene name and signed PFC for each comparison
  set.seed(100)
  pathways <- fgsea(pathways = marklist, stats = ranked_degs) # Run GSEA on DEGs for current comparison - default implementation is fgseaMultilevel w/ nPermSimple = 1000
  pathways <- select(pathways, -leadingEdge, -ES, -log2err)
  
  # Add data for states without overlap
  states_add <- unique(markers$state[markers$state %notin% pathways$pathway])
  df_add <- data.frame(pathway = states_add, pval = rep(NA, length(states_add)),
                       padj = rep(NA, length(states_add)), NES = rep(0, length(states_add)), size = rep(0, length(states_add)))
  pathways <- rbind(pathways, df_add)
  
  write.csv(pathways, paste0(output_folder, comp, "_results.csv"), row.names = FALSE)
  gsealist[[comp]] <- pathways
}

# Combined dot plot 

# Load data for all comparisons
comps <- c("nAD_NNC", "iAD_nAD", "lim_nAD", "ext_nAD")
titles <- c("nAD vs. NNC", "iAD vs. nAD", "iAD-Lim vs. nAD", "iAD-Ext vs. nAD")
names(titles) <- comps
df <- data.frame() 
for (comp in comps) {
  results <- read.csv(paste0(output_folder, comp, "_results.csv"))
  results$comp <- titles[[comp]]
  results$sig <- "Not Significant"
  results$sig[results$padj < 0.05] <- "Significant"
  results$plt_color <- paste0(results$comp, ": ", results$sig)
  df <- rbind(df, results)
}
df$plt_color <- factor(df$plt_color, levels = unique(df$plt_color)) # already ordered logically

# Order pathways by minimum NES
summary <- df %>% group_by(pathway) %>% summarize(min = min(NES))
summary <- summary %>% arrange(min)
df$pathway <- factor(df$pathway, levels = summary$pathway)

# Calculate -log10(padj) for dot size 
df$padj[is.na(df$padj)] <- 1
df$neg_logBH <- -log10(df$padj)

# Define comparison colors
colors <- c("#333333", "#303030", "#28DCF9", "#0DCDFB", "#954AFF", "#9133FF", "#FF8C3E", "#FF812C")
names(colors) <- levels(df$plt_color)

# Plot NES for all groups
plt <- ggplot() +
  geom_point(data = df, stat="identity", stroke = 1, shape = 21, fill = NA, aes(color = plt_color, x = NES, y = pathway, size = neg_logBH)) +
  scale_color_manual(values = colors) +
  labs(color = NULL) + ggtitle("GSEA Analysis of Gray Matter Microglia DEGs") +
  theme_bw() + ylab("") + xlab("NES") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_size_continuous(range = c(3, 12))

pdf(paste0(output_folder, "combined_plot.pdf"), height = 8, width = 8)
print(plt)
dev.off()


