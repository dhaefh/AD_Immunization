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
# Summary: Differential expression with MAST for amyloid-rich gray matter
#
#-----------------------------------------------

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
output_folder <- "/path/to/amyloid/mast/output/folder/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Load integrated cohort 1 Seurat object
s <- readRDS("/path/to/integrated/cohort1/object.rds")

# Define filter operator
`%notin%` <- Negate(`%in%`)

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Filter for amyloid-rich spots in gray matter for iAD and nAD, excluding vascular amyloid rich spots + neighbors
s <- subset(s, manual_annotation %notin% c("white", "meninges") & amyloid_filter == "include" & vessel_neighbor == "not_vessel" &
              amyloid_fluo > 183  & condition %in% c("iAD", "nAD"))
gc()

# Remove samples with 1 spot (Seurat error subsetting + in PrepSCTFindMarkers)
summary <- table(s$sample_id) %>% data.frame()
samples_keep <- as.character(summary$Var1[summary$Freq >= 2])
s <- subset(s, sample_id %in% samples_keep)
gc()

# Downsample spots with highest amyloid density such that no sample makes up more than 50% of a group, fold difference < 3 between groups, no more than 3000 spots per group
n_cells <- readxl::read_xlsx("/path/to/custom/downsampling/cells/per/donor.xlsx")
colnames(n_cells)[3] <- "downsampled"
cells_keep <- c()
for (sample in unique(s$sample_id)) {
  meta <- s@meta.data[s$sample_id == sample,]
  meta <- meta %>% dplyr::arrange(desc(amyloid_fluo))
  cur_cells <- n_cells$downsampled[n_cells$sample_id == sample]
  cells <- rownames(meta)[1:cur_cells]
  cells_keep <- c(cells_keep, cells)
}

# Subset for downsampled spots 
s <- subset(s, cells = cells_keep)
gc()

# Recorrect SCT data 
for (name in names(s@assays$SCT@SCTModel.list)) {
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}
s <- PrepSCTFindMarkers(s)

# Calculate CDR 
cdr <- colMeans(GetAssayData(s, assay = "SCT", layer = "data") > 0)
s$cdr <- cdr

# Standardize continuous covariates and make categorical covariates factors
s$age_centered <- scale(s$age)
s$cdr_centered <- scale(s$cdr)
s$gDNA_centered <- scale(s$gDNA_percent)
s$sex <- factor(s$sex)
s$sample_id <- factor(s$sample_id)
s$manual_annotation <- factor(s$manual_annotation)

# Define comparisons to run 
comparisons <- c("ext_vs_nAD", "lim_vs_nAD")

# Extract recorrected SCT expression data
expressionmat_full <- GetAssayData(s, assay = "SCT", layer = 'data')

# Run MAST for each comparison 
for (comparison in comparisons) {
  
  # Create subfolders
  dir.create(paste0(output_folder, comparison), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(output_folder, comparison, "/results"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(output_folder, comparison, "/plots"), showWarnings = FALSE, recursive = TRUE)
  
  # Extract names of comparison groups
  ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
  ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
  
  # If comparing iAD-Ext or iAD-Lim, set idents to condition_clearance
  if ((ident.1 %in% c("ext", "lim")) | (ident.2 %in% c("ext", "lim"))) {
    Idents(s) <- "condition_clearance"
  } else {
    Idents(s) <- "condition"
  }
  
  # Calculate fold change and percent expression using SCT data
  LFC <- FoldChange(object = s, ident.1 = ident.1, ident.2 = ident.2, fc.name = "avg_log2FC", base = 2,
                    assay = "SCT", layer = "data") %>% data.frame()
  
  # Test genes expressed in 1% of both groups and 10% of either group
  genes_keep <- row.names(LFC)[(LFC$pct.1 >= 0.01 & LFC$pct.2 >= 0.01) & (LFC$pct.1 >= 0.1 | LFC$pct.2 >= 0.1)]
  
  # Exclude contamination genes
  if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)) != 0) {
    genes_keep <- genes_keep[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)] 
  }
  
  # Filter expression matrix for genes to test
  expressionmat <- expressionmat_full[genes_keep,]
  expressionmat <- as.matrix(expressionmat)
  
  # Generate cell-level and feature-level meta data 
  fdat <- data.frame(primerid = row.names(expressionmat), row.names = row.names(expressionmat))
  cdat <- s@meta.data
  cdat$wellKey <- row.names(cdat)
  
  # Update condition variable if comparing iAD-Ext or iAD-Lim
  if ((ident.1 %in% c("ext", "lim")) | (ident.2 %in% c("ext", "lim"))) {
    cdat$condition <- cdat$condition_clearance
  } 
  
  # Create SingleCellAssay object
  sca <- FromMatrix(expressionmat, cdat, fdat, check_sanity = FALSE)
  
  # Set reference level for condition
  cond <- factor(colData(sca)$condition)
  cond <- relevel(cond, ident.1)
  colData(sca)$condition <- cond
  
  # Fit model and run LRT
  zlm <- zlm(~ condition + manual_annotation + sex + age_centered + cdr_centered + gDNA_centered + (1 | sample_id), sca, 
             ebayes = FALSE, method = "glmer")
  lrt_name <- paste0("condition", ident.2)
  summary_condition <- summary(zlm, doLRT = lrt_name) 
  
  # Extract hurdle p values
  summary_data <- summary_condition$datatable %>% data.frame()
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
  write.csv(results, paste0(output_folder, comparison, "/results/amyloid_rich.csv"))
}

