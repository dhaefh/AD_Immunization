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
#-------------------------------------------------------------------------------

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
output_folder <- "/path/to/amyloid/by/region/mast/output/folder/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Load integrated cohort 5/7/8 Seurat object
s <- readRDS("/path/to/integrated/cohort578/object.rds")

# Define variable for gDNA % 
s$gDNA_percent <- NA
s$gDNA_percent[s$sample_id == "NMA22.A1"] <- 1.4
s$gDNA_percent[s$sample_id == "NMA22.A3"] <- 1.3
s$gDNA_percent[s$sample_id == "NMA22.A4"] <- 2.5
s$gDNA_percent[s$sample_id == "NMA22.A9"] <- 1.8
s$gDNA_percent[s$sample_id == "NMA22.B1"] <- 1.0
s$gDNA_percent[s$sample_id == "NMA22.B3"] <- 0.1
s$gDNA_percent[s$sample_id == "NMA22.B4"] <- 0.2
s$gDNA_percent[s$sample_id == "NMA22.B9"] <- 0.5
s$gDNA_percent[s$sample_id == "A14.193.1"] <- 15.3
s$gDNA_percent[s$sample_id == "A14.193.3"] <- 16.3
s$gDNA_percent[s$sample_id == "A14.193.4"] <- 10.1
s$gDNA_percent[s$sample_id == "A14.193.9"] <- 1.8
s$gDNA_percent[s$sample_id == "A11.170.1"] <- 4.4
s$gDNA_percent[s$sample_id == "A11.170.3"] <- 4.9
s$gDNA_percent[s$sample_id == "A11.170.4"] <- 0.6
s$gDNA_percent[s$sample_id == "A11.170.9"] <- 0.8
print(unique(s@meta.data[,c("sample_id", "gDNA_percent")]))

# Define cortical amyloid enrichment
s$amyloid_rich <- "not_rich"
s$amyloid_rich[s$cortical_amyloid > 183 & s$vascular_amyloid == 0] <- "rich"

# Subset for cortical amyloid-rich spots in gray matter 
gray_layers <- unique(s$manual_layer[grep("gray", s$manual_layer)])
s <- subset(s, amyloid_rich == "rich" & manual_layer %in% gray_layers)
gc()

# Downsample spots with highest amyloid density

# Ensure no CAA donor makes up more than 50% of the CAA group within each region
s$caa_merged <- s$sample_id
s$caa_merged[s$caa_merged %in% c("A11.170.1", "A14.193.1")] <- "A1_new"
s$caa_merged[s$caa_merged %in% c("A11.170.3", "A14.193.3")] <- "A3_new"
s$caa_merged[s$caa_merged %in% c("A11.170.4", "A14.193.4")] <- "A4_new"
s$caa_merged[s$caa_merged %in% c("A11.170.9", "A14.193.9")] <- "A9_new"
meta <- s@meta.data
old_caa <- c()
new_caa <- c()
for (region in c(1, 3, 4, 9)) {
  
  # If old CAA donor has more spots, downsample to total new controls, otherwise ensure no new donors not make up more than 50% of the total CAA pool
  if (sum(meta$caa_merged == paste0("NMA22.A", region)) > sum(meta$caa_merged == paste0("A", region, "_new"))) {
    cur_meta <- meta[meta$caa_merged == paste0("NMA22.A", region),]
    cur_meta <- cur_meta %>% dplyr::arrange(desc(cortical_amyloid))
    cur_downsample <- sum(meta$caa_merged == paste0("A", region, "_new"))
    cells <- rownames(cur_meta)[1:cur_downsample]
    old_caa <- c(old_caa, cells)
    new_caa <- c(new_caa, rownames(meta)[meta$caa_merged == paste0("A", region, "_new")])
  } else {
    
    # Calculate spots per new CAA donor
    total_pool <- sum(meta$caa_merged %in% c(paste0("A", region, "_new"), paste0("NMA22.A", region))) 
    caa_summary <- meta[meta$caa_merged == paste0("A", region, "_new"),] %>% group_by(sample_id) %>% summarize(count = n()) 
    
    # Take ceiling if downsampling to be conservative
    if (max(caa_summary$count) <= round(ceiling(total_pool/2))) {
      new_caa <- c(new_caa, rownames(meta)[meta$caa_merged == paste0("A", region, "_new")])
      old_caa <- c(old_caa, rownames(meta)[meta$caa_merged == paste0("NMA22.A", region)])
    } else {
      sample <- caa_summary$sample_id[caa_summary$count > round(ceiling(total_pool/2))]
      n_cells <- sum(meta$caa_merged == paste0("NMA22.A", region) | (meta$caa_merged == paste0("A", region, "_new") & meta$sample_id != sample)) # Calculate remaining spots in the CAA pool
      
      cur_meta <- meta[meta$sample_id == sample,]
      cur_meta <- cur_meta %>% dplyr::arrange(desc(cortical_amyloid))
      cells <- rownames(cur_meta)[1:n_cells]
      
      new_caa <- c(new_caa, rownames(meta)[rownames(meta) %in% cells | (meta$caa_merged == paste0("A", region, "_new") & meta$sample_id != sample)])
      old_caa <- c(old_caa, rownames(meta)[meta$caa_merged == paste0("NMA22.A", region)])
    }
  }
}

# Filter for downsampled controls
meta <- meta[rownames(meta) %in% c(old_caa, new_caa) | str_detect(meta$sample_id, "B"),]

# Define variable for total CAA per region
meta$caa_merged_total <- meta$caa_merged
meta$caa_merged_total[grep("A1", meta$caa_merged_total)] <- "A1"
meta$caa_merged_total[grep("A3", meta$caa_merged_total)] <- "A3"
meta$caa_merged_total[grep("A4", meta$caa_merged_total)] <- "A4"
meta$caa_merged_total[grep("A9", meta$caa_merged_total)] <- "A9"

# Within each region, ensure fold difference < 3 between groups and each group has no more than 3000 spots
total_cells_keep <- c()
for (region in c(1, 3, 4, 9)) {
  
  # Identify group with more spots
  if (sum(meta$caa_merged_total == paste0("A", region)) > sum(meta$sample_id == paste0("NMA22.B", region))) {
    max_group <- "CAA"
    max <- sum(meta$caa_merged_total == paste0("A", region))
    min <- sum(meta$sample_id == paste0("NMA22.B", region))
  } else {
    max_group <- "LCMB"
    max <- sum(meta$sample_id == paste0("NMA22.B", region))
    min <- sum(meta$caa_merged_total == paste0("A", region))
  }
  
  # Calculate fold difference
  fold <- max/min
  
  # If fold difference > 3 and CAA is the larger group, downsample CAA to 3x LCMB
  if (fold > 3 & max_group == "CAA") {
    target <- min*3
    
    # Adjust if greater than 3000
    if (target > 3000) {
      target <- 3000
    }
    
    # Calculate z value 
    z <- max - target
    
    # Calculate target spots per donor
    target_per_donor <- round(target/length(unique(meta$sample_id[meta$caa_merged_total == paste0("A", region)])))
    
    # Calculate deviation from target per donor
    caa_summary <- meta[meta$caa_merged_total == paste0("A", region),] %>% group_by(sample_id) %>% summarize(count = n())
    caa_summary$deviation <- caa_summary$count - target_per_donor
    
    # Arrange by decreasing deviation
    caa_summary <- caa_summary %>% arrange(desc(deviation))
    
    # Filter for donors with positive deviation
    positive_summary <- caa_summary[caa_summary$deviation > 0,]
    
    # If z is large enough to reach target per donor, downsample all to target per donor
    if (z >= sum(positive_summary$deviation)) {
      positive_summary$n_sample <- target_per_donor
    } 
    
    # Else if z is large enough to reach min positive deviation, downsample to min positive deviation 
    else if (z >= sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation))) {
      
      positive_summary$n_sample <- positive_summary$count[positive_summary$deviation == min(positive_summary$deviation)]
      
      # If possible, downsample positive donors equally using remaining z
      remaining_z <- z - sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation)) # Subtract spots removed in previous step from total z
      if (remaining_z > 0) {
        positive_summary$n_sample <- positive_summary$n_sample - round(remaining_z/nrow(positive_summary))
      }
    }
    
    # Else if z is at least the difference between top 2 positive deviations, equalize the top 2, then downsample top 2 equally 
    else if (z >= positive_summary$deviation[1] - positive_summary$deviation[2]) {
      
      positive_summary$n_sample <- positive_summary$count # Initialize with actual count
      positive_summary$n_sample[1] <- positive_summary$n_sample[2] # Downsample highest to second highest
      
      # If possible, downsample top 2 donors equally using remaining z
      remaining_z <- z - (positive_summary$deviation[1] - positive_summary$deviation[2])
      if (remaining_z > 0) {
        positive_summary$n_sample[1:2] <- positive_summary$n_sample[1:2] - round(remaining_z/2)
      }
    } 
    
    # Otherwise just downsample largest donor
    else if (z < positive_summary$deviation[1] - positive_summary$deviation[2]) {
      positive_summary$n_sample <- positive_summary$count
      positive_summary$n_sample[1] <- positive_summary$n_sample[1] - z
    } else {
      print("Case not covered")
    }
    
    # Downsample the positive donors, keep all spots for the others
    cells_keep <- c()
    for (sample in unique(meta$sample_id[meta$caa_merged_total == paste0("A", region)])) {
      if (sample %in% positive_summary$sample_id) {
        
        cur_downsample <- positive_summary$n_sample[positive_summary$sample_id == sample]
        cur_meta <- meta[meta$sample_id == sample,]
        cur_meta <- cur_meta %>% dplyr::arrange(desc(cortical_amyloid))
        cells <- rownames(cur_meta)[1:cur_downsample]
        
      } else {
        cells <- rownames(meta)[meta$sample_id == sample]
      }
      cells_keep <- c(cells_keep, cells)
    }
    
    # Keep all LCMB + downsampled CAA 
    cells_keep <- rownames(meta)[rownames(meta) %in% cells_keep | meta$sample_id == paste0("NMA22.B", region)]
  }
  # If fold difference > 3 and LCMB is the larger group, downsample LCMB to 3x CAA
  else if (fold > 3 & max_group == "LCMB") { 
    target <- min*3
    
    # Adjust if greater than 3000
    if (target > 3000) {
      target <- 3000
    }
    cur_meta <- meta[meta$sample_id == paste0("NMA22.B", region),]
    cur_meta <- cur_meta %>% dplyr::arrange(desc(cortical_amyloid))
    cells_keep <- rownames(cur_meta)[1:target]
    cells_keep <- rownames(meta)[rownames(meta) %in% cells_keep | meta$caa_merged_total == paste0("A", region)] 
  } 
  
  # Downsample CAA if fold difference <= 3 but CAA has > 3000 spots (this case only applies for CAA in this cohort)
  else if (sum(meta$caa_merged_total == paste0("A", region)) > 3000) {
    target <- 3000
    
    # Calculate z value for CAA 
    z <- sum(meta$caa_merged_total == paste0("A", region)) - target
    
    # Calculate target spots per donor
    target_per_donor <- round(target/length(unique(meta$sample_id[meta$caa_merged_total == paste0("A", region)])))
    
    # Calculate deviation from target per donor
    caa_summary <- meta[meta$caa_merged_total == paste0("A", region),] %>% group_by(sample_id) %>% summarize(count = n())
    caa_summary$deviation <- caa_summary$count - target_per_donor
    
    # Arrange by decreasing deviation
    caa_summary <- caa_summary %>% arrange(desc(deviation))
    
    # Filter for donors with positive deviation
    positive_summary <- caa_summary[caa_summary$deviation > 0,]
    
    # If z is large enough to reach target per donor, downsample all to target per donor
    if (z >= sum(positive_summary$deviation)) {
      positive_summary$n_sample <- target_per_donor
    } 
    
    # Else if z is large enough to reach min positive deviation, downsample to min positive deviation 
    else if (z >= sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation))) {
      
      positive_summary$n_sample <- positive_summary$count[positive_summary$deviation == min(positive_summary$deviation)]
      
      # If possible, downsample positive donors equally using remaining z
      remaining_z <- z - sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation)) # Subtract spots removed in previous step from total z
      if (remaining_z > 0) {
        positive_summary$n_sample <- positive_summary$n_sample - round(remaining_z/nrow(positive_summary))
      }
    }
    
    # Else if z is at least the difference between top 2 positive deviations, equalize the top 2, then downsample top 2 equally 
    else if (z >= positive_summary$deviation[1] - positive_summary$deviation[2]) {
      
      positive_summary$n_sample <- positive_summary$count # Initialize with actual count
      positive_summary$n_sample[1] <- positive_summary$n_sample[2] # Downsample highest to second highest
      
      # If possible, downsample top 2 donors equally using remaining z
      remaining_z <- z - (positive_summary$deviation[1] - positive_summary$deviation[2])
      if (remaining_z > 0) {
        positive_summary$n_sample[1:2] <- positive_summary$n_sample[1:2] - round(remaining_z/2)
      }
    } 
    
    # Otherwise just downsample largest donor
    else if (z < positive_summary$deviation[1] - positive_summary$deviation[2]) {
      positive_summary$n_sample <- positive_summary$count
      positive_summary$n_sample[1] <- positive_summary$n_sample[1] - z
    } else {
      print("Case not covered")
    }
    
    # Downsample the positive donors, keep all spots for the others
    cells_keep <- c()
    for (sample in unique(meta$sample_id[meta$caa_merged_total == paste0("A", region)])) {
      if (sample %in% positive_summary$sample_id) {
        cur_downsample <- positive_summary$n_sample[positive_summary$sample_id == sample]
        cur_meta <- meta[meta$sample_id == sample,]
        cur_meta <- cur_meta %>% dplyr::arrange(desc(cortical_amyloid))
        cells <- rownames(cur_meta)[1:cur_downsample]
      } else {
        cells <- rownames(meta)[meta$sample_id == sample]
      }
      cells_keep <- c(cells_keep, cells)
    }
    
    # Keep all LCMB + downsampled CAA 
    cells_keep <- rownames(meta)[rownames(meta) %in% cells_keep | meta$sample_id == paste0("NMA22.B", region)]
  } else { 
    cells_keep <- rownames(meta)[meta$caa_merged_total == paste0("A", region) | meta$sample_id == paste0("NMA22.B", region)] 
  }
  total_cells_keep <- c(total_cells_keep, cells_keep)
}

# Filter for final downsampled spots
meta <- meta[total_cells_keep,]


# Generate table of spots per donor, ordered by region 
df <- data.frame(table(meta$sample_id))
colnames(df) <- c("sample_id", "downsampled")
rownames(df) <- df$sample_id
df <- df[c(paste0(c("A11.170.", "A14.193."), "1"), "NMA22.A1", "NMA22.B1",
           paste0(c("A11.170.", "A14.193."), "3"), "NMA22.A3", "NMA22.B3",
           paste0(c("A11.170.", "A14.193."), "4"), "NMA22.A4", "NMA22.B4",
           paste0(c("A11.170.", "A14.193."), "9"), "NMA22.A9", "NMA22.B9"),]
df$sample_id <- NULL

df2 <- data.frame(table(s$sample_id))
colnames(df2) <- c("sample_id", "raw")
rownames(df2) <- df2$sample_id
df2 <- df2[rownames(df),]
df$raw <- df2$raw
s
write.csv(df[,c("raw", "downsampled")], paste0(output_folder, "downsampling_final.csv"))

# Subset for downsampled spots
s <- subset(s, cells = rownames(meta))
gc()

# Create group variable
s$group_de <- s$sample_id
s$group_de[grep("\\.A1|193.1|170.1", s$group_de)] <- "A1"
s$group_de[grep("\\.A3|193.3|170.3", s$group_de)] <- "A3"
s$group_de[grep("\\.A4|193.4|170.4", s$group_de)] <- "A4"
s$group_de[grep("\\.A9|193.9|170.9", s$group_de)] <- "A9"
s$group_de[grep("\\.B1", s$group_de)] <- "B1"
s$group_de[grep("\\.B3", s$group_de)] <- "B3"
s$group_de[grep("\\.B4", s$group_de)] <- "B4"
s$group_de[grep("\\.B9", s$group_de)] <- "B9"

# Standardize continuous covariates and make categorical covariates factors
s$group_de <- factor(s$group_de)
s$sample_id <- factor(s$sample_id)
s$manual_layer <- factor(s$manual_layer)
s$gDNA_centered <- scale(s$gDNA_percent)

# Recorrect SCT data 
for (name in names(s@assays$SCT@SCTModel.list)) {
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}
s <- PrepSCTFindMarkers(s)

# Calculate CDR 
cdr <- colMeans(GetAssayData(s, assay = "SCT", layer = "data") > 0)
s$cdr <- cdr

# Standardize CDR
s$cdr_centered <- scale(s$cdr)

# Define DE comparisons to run 
comparisons <- c("B1_vs_A1", "B3_vs_A3", "B4_vs_A4", "B9_vs_A9")

# Extract recorrected SCT expression data
expressionmat_full <- GetAssayData(s, assay = "SCT", layer = 'data')

# Set idents 
Idents(s) <- "group_de"

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Run MAST for each comparison
for (comparison in comparisons) {
  
  # Create output folders
  dir.create(paste0(output_folder, comparison), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(output_folder, comparison, "/results"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(output_folder, comparison, "/plots"), showWarnings = FALSE, recursive = TRUE)
  
  # Extract names of comparison groups
  ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
  ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
  
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
  
  # Create SingleCellAssay object
  sca <- FromMatrix(expressionmat, cdat, fdat, check_sanity = FALSE)
  
  # Set reference level for comparison group
  group <- factor(colData(sca)$group_de)
  group <- relevel(group, ident.1) 
  colData(sca)$group_de <- group
  
  # Fit model and run LRT
  zlm_group <- zlm(~ group_de + manual_layer + gDNA_centered + cdr_centered + (1 | sample_id), sca, ebayes = FALSE, method = "glmer")
  lrt_name <- paste0("group_de", ident.2)
  summary <- summary(zlm_group, doLRT = lrt_name) # FindMarkers: doLRT = "conditionGroup2"
  
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
  write.csv(results, paste0(output_folder, comparison, "/results/amyloid_rich.csv"))
}

