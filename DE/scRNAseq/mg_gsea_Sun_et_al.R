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
# Summary: Microglia state enrichment analysis with fgsea, using human microglia states from Sun et al. 
#
#-----------------------------------------------
# Initialization

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
output_folder <- "/path/to/Sun_et_al/gsea/output/folder/"

# Load microglia state markers
markers <- read_excel("/path/to/Sun_et_al/states.xlsx")

# Filter operator 
`%notin%` <- Negate(`%in%`)

# Update column names 
markers <- markers %>% rename(state = microgliaState)

# Group markers by microglia state
markers <- group_by(markers, state)

# Get list of state markers
marklist <- lapply(split(select(markers, gene), markers$state), deframe)

# Run GSEA for all regions

# Load results and calculate signed PFC, ensuring values are finite 
results <- read.csv("/path/to/all/regions/celltype/mast/output/results/Mg.csv", row.names = 1)
results$PFC <- -log10(results$BH)*results$avg_log2FC
print(sum(results$PFC == Inf))

# Run GSEA
options(mc.cores = 1)
ranked_degs <- deframe(select(results, gene, PFC)) # Get gene name and signed PFC 
set.seed(100)
pathways <- fgsea(pathways = marklist, stats = ranked_degs) # Run GSEA on DEGs - default implementation is fgseaMultilevel w/ nPermSimple = 1000
pathways <- select(pathways, -leadingEdge, -ES, -log2err)

# Add data for states without overlap
states_add <- unique(markers$state[markers$state %notin% pathways$pathway])
df_add <- data.frame(pathway = states_add, pval = rep(NA, length(states_add)),
                     padj = rep(NA, length(states_add)), NES = rep(0, length(states_add)), size = rep(0, length(states_add)))
pathways <- rbind(pathways, df_add)

# Save results
write.csv(pathways, paste0(output_folder, "all_regions_results.csv"), row.names = FALSE)

# Run GSEA per region 
de_folder <- "/path/to/by/region/celltype/mast/output/"
for (comp in c("B1_vs_A1", "B3_vs_A3", "B4_vs_A4", "B9_vs_A9")) {
  
  # Load results and calculate signed PFC for DEGs, ensuring values are finite
  results <- read.csv(paste0(de_folder, comp, "/results/Mg.csv"), row.names = 1)
  results$PFC <- -log10(results$BH)*results$avg_log2FC
  print(sum(results$PFC == Inf))
  
  # Define region for comparison  
  if (str_detect(comp, "1")) {
    region <- "FCX"
  } else if (str_detect(comp, "3")) {
    region <- "TCX"
  } else if (str_detect(comp, "4")) {
    region <- "PCX"
  } else {
    region <- "HIPP"
  }
  
  # Run GSEA
  options(mc.cores = 1)
  ranked_degs <- deframe(select(results, gene, PFC)) # Get gene name and signed PFC 
  set.seed(100)
  pathways <- fgsea(pathways = marklist, stats = ranked_degs) # Run GSEA on DEGs - default implementation is fgseaMultilevel w/ nPermSimple = 1000
  pathways <- select(pathways, -leadingEdge, -ES, -log2err)
  
  # Add data for states without overlap
  states_add <- unique(markers$state[markers$state %notin% pathways$pathway])
  df_add <- data.frame(pathway = states_add, pval = rep(NA, length(states_add)),
                       padj = rep(NA, length(states_add)), NES = rep(0, length(states_add)), size = rep(0, length(states_add)))
  pathways <- rbind(pathways, df_add)
  
  # Save results
  print(sum(is.na(pathways)))
  write.csv(pathways, paste0(output_folder, region, "_results.csv"), row.names = FALSE)
  
  # Generate bar plot
  plt_data <- pathways 
  plt_data$significance[is.na(plt_data$NES)] <- "Could not Compute P Value" # Label pathways where p value/NES could not be calculated
  plt_data$significance[is.na(plt_data$padj) & !is.na(plt_data$NES)] <- "No Gene Overlap"
  plt_data$NES[is.na(plt_data$NES)] <- 0 # Set NES to zero if NA
  plt_data$significance[plt_data$padj < 0.05] <- "Significant (padj < 0.05)"
  plt_data$significance[plt_data$padj >= 0.05] <- "Not Significant"
  plt <- ggplot(plt_data, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=significance)) + coord_flip() +
    theme_cowplot() + ggtitle(paste0("LCMB vs. CAA: ", region)) +
    scale_fill_manual(values = c("Not Significant" = "#C77CFF", "Significant (padj < 0.05)" = "#FDB863", "No Gene Overlap" = "black",
                                 "Could not Compute P Value" = "gray")) + labs(x = NULL, fill = NULL)
  pdf(paste0(output_folder, region, "_plot.pdf"), height = 8, width = 8) 
  print(plt)
  dev.off()
}


# Combined dot plot 

# Load data for all comparisons
comps <- c("FCX", "TCX", "PCX", "HIPP", "all_regions")
titles <- c("FCX", "TCX", "PCX", "HIPP", "All Regions")
names(titles) <- comps
df <- data.frame() 
for (comp in comps) {
  results <- read.csv(paste0(output_folder, comp, "_results.csv"))
  results$comp <- titles[[comp]]
  results$sig <- "Not Significant"
  results$sig[is.na(results$NES)] <- "Could not Compute P Value" # Label pathways where p value/NES could not be calculated
  results$NES[is.na(results$NES)] <- 0 # Set NES to zero if NA
  results$sig[results$padj < 0.05] <- "Significant"
  results$plt_color <- paste0(results$comp, ": ", results$sig)
  df <- rbind(df, results)
}
df$plt_color <- factor(df$plt_color, levels = unique(df$plt_color)) 

# Order pathways by minimum NES
summary <- df %>% group_by(pathway) %>% summarize(min = min(NES))
summary <- summary %>% arrange(min)
df$pathway <- factor(df$pathway, levels = summary$pathway)

# Calculate -log10(padj) for dot size 
df$padj[is.na(df$padj)] <- 1
df$neg_logBH <- -log10(df$padj)

# Define comparison colors
colors <- c("#D682FF", "#C96CFF", "#68E7C4", "#EC2E8B", "#EA1F8A", "gray", "#FFBF5F", "#FFB14A", "#333333", "#303030")
names(colors) <- levels(df$plt_color)

# Plot NES for all groups
plt <- ggplot() +
  geom_point(data = df, stat="identity", stroke = 1, shape = 21, fill = NA, aes(color = plt_color, x = NES, y = pathway, size = neg_logBH)) +
  scale_color_manual(values = colors) +
  labs(color = NULL) + ggtitle("GSEA Analysis of Microglia DEGs: LCMB vs. CAA") +
  theme_bw() + ylab("") + xlab("NES") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_size_continuous(range = c(3, 12))

pdf(paste0(output_folder, "combined_plot.pdf"), height = 8, width = 8)
print(plt)
dev.off()








