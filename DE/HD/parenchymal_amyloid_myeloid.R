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
# Summary: Differential expression with MAST for parenchymal amyloid associated myeloid cells
#
#-----------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("DESeq2")
  library("ggrepel")
  library("UpSetR")
  library("randomcoloR")
  library("ComplexUpset")
  library("MAST")
})

# Define output folder
hd_folder <- "/path/to/general/HD/output/folder/"
output_folder <- paste0(hd_folder, "differential_expression/mast/")
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_folder, "between_groups"), showWarnings = FALSE, recursive = TRUE)

# Load integrated HD object
s <- readRDS(paste0(hd_folder, "integration/data/s_integrated.rds"))

# Define parenchymal amyloid distance threshold in pixels (threshold = 20 micron, there are ~0.50292 micron per pixel)
amyloid_thresh <- 20/0.50292

# Identify parenchymal amyloid-associated myeloid cells in hippocampus 
s$keep <- FALSE
s$keep[s$hippocampus == TRUE & s$parenchymal_amyloid_distance < amyloid_thresh & s$harmony_snn_res.0.2 == 5] <- TRUE
s <- subset(s, keep == TRUE)
table(s$sample_id)

# Define variable for DE group 
s$group_de <- ifelse(s$sample_id == "B9", "LCMB", "nAD")
unique(s@meta.data[,c("sample_id", "group_de")])
Idents(s) <- "group_de"

# Recorrect SCT data
DefaultAssay(s) <- "SCT"
s <- PrepSCTFindMarkers(s)

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5)

# Calculate CDR
cdr <- colMeans(GetAssayData(s, assay = "SCT", layer = "data") > 0)
s$cdr <- cdr

# Standardize CDR 
s$cdr_centered <- scale(s$cdr)

# Make sample ID a factor
s$sample_id <- factor(s$sample_id)

# Extract recorrected SCT expression data
expressionmat_full <- GetAssayData(s, assay = "SCT", layer = 'data')

# Define comparison to run and extract names of comparison groups
comparison <- "LCMB_vs_nAD"
ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]

# Calculate fold change and percent expression using SCT data
LFC <- FoldChange(object = s, ident.1 = ident.1, ident.2 = ident.2, fc.name = "avg_log2FC", base = 2,
                  assay = "SCT", layer = "data") %>% data.frame()

# Test genes expressed in 1% of either group
genes_keep <- row.names(LFC)[LFC$pct.1 >= 0.01 | LFC$pct.2 >= 0.01]

# Exclude contamination genes
if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)) != 0) {
  genes_keep <- genes_keep[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)] 
}

# Filter expression matrix for genes to test
expressionmat <- expressionmat_full[genes_keep,]
expressionmat <- as.matrix(expressionmat)

# Generate cell-level and feature-level meta data 
cdat <- s@meta.data
cdat$wellKey <- row.names(cdat)
fdat <- data.frame(primerid = row.names(expressionmat), row.names = row.names(expressionmat))

# Create SingleCellAssay object
sca <- FromMatrix(expressionmat, cdat, fdat, check_sanity = FALSE)

# Set reference level for DE group
group <- factor(colData(sca)$group_de)
group <- relevel(group, ident.1) 
colData(sca)$group_de <- group

# Fit model and run LRT
options(mc.cores = 1)
zlm_group <- zlm(~ group_de + cdr_centered, sca)
lrt_name <- paste0("group_de", ident.2)
summary <- summary(zlm_group, doLRT = lrt_name) 

# Extract hurdle p values
summary_data <- summary$datatable %>% data.frame()
p_val <- summary_data[summary_data[, "component"] == "H", 4]
genes.return <- summary_data[summary_data[, "component"] == "H", 1]

# Compile results
results <- data.frame(p_val = p_val, gene = genes.return, row.names = genes.return)
results$BH <- p.adjust(results$p_val, method = "BH")
LFC_use <- LFC[match(results$gene, row.names(LFC)),]
print(sum(row.names(LFC_use) != results$gene))
results$avg_log2FC <- LFC_use$avg_log2FC

# Define DEGs
results$DE <- "Not DE"
results$DE[results$avg_log2FC > fc_thresh & results$BH < p_thresh] <- "Upregulated"
results$DE[results$avg_log2FC < -fc_thresh & results$BH < p_thresh] <- "Downregulated"
results$DE_gene <- NA
results$DE_gene[results$DE != "Not DE"] <- results$gene[results$DE != "Not DE"]

# Save results
write.csv(results, paste0(output_folder, "between_groups/results.csv"))

