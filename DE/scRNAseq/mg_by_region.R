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
# Summary: Differential expression with MAST for microglia 
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
  library("MAST")
})

# Define output folder
output_folder <- "/path/to/by/region/celltype/mast/output/"

# Load integrated scRNAseq Seurat object
s <- readRDS("/path/to/integrated/scRNAseq/object.rds")

# Subset for immune subtypes
immune <- readRDS("/path/to/integrated/scRNAseq/immune/object.rds")
s <- subset(s, cells = rownames(immune@meta.data))
immune <- NULL
gc()

# Define group variable
s@meta.data$group_de <- NA
s$group_de[grep("^102", s$sample_merged)] <- "AN1792"
s$group_de[grep("\\.B", s$sample_merged)] <- "LCMB"
s$group_de[grep("\\.A", s$sample_merged)] <- "CAA"

# Remove AN1792 samples
s <- subset(s, group_de != "AN1792")
gc()

# Format cell type names
s$cell_type_de <- s$merged_celltype_final
s$cell_type_de[s$cell_type_de == "Infl. Endo"] <- "Infl_Endo"

# Define variable combining new CAA controls within each region
s$caa_merged <- s$sample_merged
s$caa_merged[s$caa_merged %in% c("A11.170.A1", "A14.193.A1")] <- "A1_new"
s$caa_merged[s$caa_merged %in% c("A11.170.A3", "A14.193.A3")] <- "A3_new"
s$caa_merged[s$caa_merged %in% c("A11.170.A4", "A14.193.A4")] <- "A4_new"
s$caa_merged[s$caa_merged %in% c("A11.170.A9", "A14.193.A9")] <- "A9_new"

# Downsample CAA and LCMB per cell type per region
full_meta <- data.frame()
for (type in unique(s$cell_type_de)) {
  
  # Extract meta data for current cell type 
  print(type)
  meta <- s@meta.data[s$cell_type_de == type,]
  
  # Skip if either group is not present 
  if (length(unique(meta$group_de)) != 2) {
    next
  }
  
  # Ensure no CAA donor makes up more than 50% of the CAA group within each region
  old_caa <- c()
  new_caa <- c()
  for (region in c("A1", "A3", "A4", "A9")) {
    
    # Skip if no CAA in region 
    if (sum(meta$caa_merged %in% c(paste0("NMA22.300.", region), paste0(region, "_new"))) == 0) {
      next
    }
    
    # If old CAA donor has more cells, downsample to total new controls, otherwise ensure no new donors make up more than 50% of the total CAA pool
    if (sum(meta$caa_merged == paste0("NMA22.300.", region)) > sum(meta$caa_merged == paste0(region, "_new"))) {
      set.seed(100)
      cells <- sample(rownames(meta)[meta$caa_merged == paste0("NMA22.300.", region)], sum(meta$caa_merged == paste0(region, "_new")), replace = FALSE)
      old_caa <- c(old_caa, cells)
      new_caa <- c(new_caa, rownames(meta)[meta$caa_merged == paste0(region, "_new")])
    } else {
      
      # Calculate cells per new CAA donor
      total_pool <- sum(meta$caa_merged %in% c(paste0(region, "_new"), paste0("NMA22.300.", region))) 
      caa_summary <- meta[meta$caa_merged == paste0(region, "_new"),] %>% group_by(sample_merged) %>% summarize(count = n()) 
      
      # Downsample if one donor makes up more than 50% (rounding up to exclude border cases) of the pool
      if (max(caa_summary$count) <= round(ceiling(total_pool/2))) {
        new_caa <- c(new_caa, rownames(meta)[meta$caa_merged == paste0(region, "_new")])
        old_caa <- c(old_caa, rownames(meta)[meta$caa_merged == paste0("NMA22.300.", region)])
      } else {
        sample <- caa_summary$sample_merged[caa_summary$count > round(ceiling(total_pool/2))]
        n_cells <- sum(meta$caa_merged == paste0("NMA22.300.", region) | (meta$caa_merged == paste0(region, "_new") & meta$sample_merged != sample)) 
        set.seed(100)
        cells <- sample(rownames(meta)[meta$sample_merged == sample], n_cells, replace = FALSE)
        new_caa <- c(new_caa, rownames(meta)[rownames(meta) %in% cells | (meta$caa_merged == paste0(region, "_new") & meta$sample_merged != sample)])
        old_caa <- c(old_caa, rownames(meta)[meta$caa_merged == paste0("NMA22.300.", region)])
      }
    }
  }
  
  # Filter for downsampled controls
  meta <- meta[rownames(meta) %in% c(old_caa, new_caa) | meta$group_de == "LCMB",]
  
  # Define variable for total CAA per region
  meta$caa_merged_total <- meta$caa_merged
  meta$caa_merged_total[grep("A1", meta$caa_merged_total)] <- "A1"
  meta$caa_merged_total[grep("A3", meta$caa_merged_total)] <- "A3"
  meta$caa_merged_total[grep("A4", meta$caa_merged_total)] <- "A4"
  meta$caa_merged_total[grep("A9", meta$caa_merged_total)] <- "A9"
  
  # Within each region, ensure fold difference < 3 between groups and each group has no more than 3000 cells
  total_cells_keep <- c()
  for (region in c(1, 3, 4, 9)) {
    
    # Skip if either group is missing
    if (sum(meta$caa_merged_total == paste0("A", region)) == 0 | sum(meta$sample_merged == paste0("NMA22.205.B", region)) == 0) {
      next
    }
    
    # If both groups have > 3000, set target max of both to 3000 and downsample accordingly
    if (sum(meta$caa_merged_total == paste0("A", region)) > 3000 & sum(meta$sample_merged == paste0("NMA22.205.B", region)) > 3000) {
      
      # Set uniform target
      target <- 3000
      
      # Calculate cells to remove for CAA
      z <- sum(meta$caa_merged_total == paste0("A", region)) - target
      
      # Calculate target cells per donor
      target_per_donor <- round(target/length(unique(meta$sample_merged[meta$caa_merged_total == paste0("A", region)])))
      
      # Calculate deviation from target per donor
      caa_summary <- meta[meta$caa_merged_total == paste0("A", region),] %>% group_by(sample_merged) %>% summarize(count = n())
      caa_summary$deviation <- caa_summary$count - target_per_donor
      
      # Arrange by decreasing deviation
      caa_summary <- caa_summary %>% arrange(desc(deviation))
      
      # Filter for donors with positive deviation
      positive_summary <- caa_summary[caa_summary$deviation > 0,]
      
      # If z is large enough to reach target per donor, downsample all to target per donor
      if (z >= sum(positive_summary$deviation)) {
        positive_summary$n_sample <- target_per_donor
      } 
      
      # Else if z is large enough to reach min positive deviation, downsample to min positive deviation, then downsample equally
      else if (z >= sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation))) {
        positive_summary$n_sample <- positive_summary$count[positive_summary$deviation == min(positive_summary$deviation)]
        remaining_z <- z - sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation)) 
        if (remaining_z > 0) {
          positive_summary$n_sample <- positive_summary$n_sample - round(remaining_z/nrow(positive_summary))
        }
      }
      
      # Else if z is at least the difference between top 2 positive deviations, equalize the top 2, then downsample top 2 equally 
      else if (z >= positive_summary$deviation[1] - positive_summary$deviation[2]) {
        positive_summary$n_sample <- positive_summary$count 
        positive_summary$n_sample[1] <- positive_summary$n_sample[2]
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
      
      # Downsample the positive donors, keep all cells for the others
      cells_keep <- c()
      for (sample in unique(meta$sample_merged[meta$caa_merged_total == paste0("A", region)])) {
        if (sample %in% positive_summary$sample_merged) {
          set.seed(100)
          cells <- sample(rownames(meta)[meta$sample_merged == sample], positive_summary$n_sample[positive_summary$sample_merged == sample], replace = FALSE)
        } else {
          cells <- rownames(meta)[meta$sample_merged == sample]
        }
        cells_keep <- c(cells_keep, cells)
      }
      
      # Downsample LCMB to target
      set.seed(100)
      cells <- sample(rownames(meta)[meta$sample_merged == paste0("NMA22.205.B", region)], target, replace = FALSE)
      cells_keep <- c(cells_keep, cells)
    } else {
      
      # Identify larger group
      if (sum(meta$caa_merged_total == paste0("A", region)) > sum(meta$sample_merged == paste0("NMA22.205.B", region))) {
        max_group <- "CAA"
        max <- sum(meta$caa_merged_total == paste0("A", region))
        min <- sum(meta$sample_merged == paste0("NMA22.205.B", region))
      } else {
        max_group <- "LCMB"
        max <- sum(meta$sample_merged == paste0("NMA22.205.B", region))
        min <- sum(meta$caa_merged_total == paste0("A", region))
      }
      
      # Calculate fold difference
      fold <- max/min
      
      # If fold difference > 3 and CAA is the larger group, downsample CAA to 3x LCMB
      if (fold > 3 & max_group == "CAA") {
        target <- min*3
        
        # Adjust target if greater than 3000
        if (target > 3000) {
          target <- 3000
        }
        
        # Calculate cells to remove 
        z <- max - target
        
        # Calculate target cells per donor
        target_per_donor <- round(target/length(unique(meta$sample_merged[meta$caa_merged_total == paste0("A", region)])))
        
        # Calculate deviation from target per donor
        caa_summary <- meta[meta$caa_merged_total == paste0("A", region),] %>% group_by(sample_merged) %>% summarize(count = n())
        caa_summary$deviation <- caa_summary$count - target_per_donor
        
        # Arrange by decreasing deviation
        caa_summary <- caa_summary %>% arrange(desc(deviation))
        
        # Filter for donors with positive deviation
        positive_summary <- caa_summary[caa_summary$deviation > 0,]
        
        # If z is large enough to reach target per donor, downsample all to target per donor
        if (z >= sum(positive_summary$deviation)) {
          positive_summary$n_sample <- target_per_donor
        } 
        
        # Else if z is large enough to reach min positive deviation, downsample to min positive deviation, then downsample equally
        else if (z >= sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation))) {
          positive_summary$n_sample <- positive_summary$count[positive_summary$deviation == min(positive_summary$deviation)]
          remaining_z <- z - sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation)) 
          if (remaining_z > 0) {
            positive_summary$n_sample <- positive_summary$n_sample - round(remaining_z/nrow(positive_summary))
          }
        }
        
        # Else if z is at least the difference between top 2 positive deviations, equalize the top 2, then downsample top 2 equally 
        else if (z >= positive_summary$deviation[1] - positive_summary$deviation[2]) {
          positive_summary$n_sample <- positive_summary$count 
          positive_summary$n_sample[1] <- positive_summary$n_sample[2] 
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
        
        # Downsample the positive donors, keep all cells for the others
        cells_keep <- c()
        for (sample in unique(meta$sample_merged[meta$caa_merged_total == paste0("A", region)])) {
          if (sample %in% positive_summary$sample_merged) {
            set.seed(100)
            cells <- sample(rownames(meta)[meta$sample_merged == sample], positive_summary$n_sample[positive_summary$sample_merged == sample], replace = FALSE)
          } else {
            cells <- rownames(meta)[meta$sample_merged == sample]
          }
          cells_keep <- c(cells_keep, cells)
        }
        
        # Keep all LCMB + downsampled CAA 
        cells_keep <- rownames(meta)[rownames(meta) %in% cells_keep | meta$sample_merged == paste0("NMA22.205.B", region)]
      }
      # If fold difference > 3 and LCMB is the larger group, downsample LCMB to 3x CAA
      else if (fold > 3 & max_group == "LCMB") { 
        target <- min*3
        
        # Adjust target if greater than 3000
        if (target > 3000) {
          target <- 3000
        }
        set.seed(100)
        cells_keep <- sample(rownames(meta)[meta$sample_merged == paste0("NMA22.205.B", region)], target, replace = FALSE)
        cells_keep <- rownames(meta)[rownames(meta) %in% cells_keep | meta$caa_merged_total == paste0("A", region)] 
      } 
      
      # Downsample CAA if fold difference <= 3 but CAA has > 3000 cells
      else if (sum(meta$caa_merged_total == paste0("A", region)) > 3000) {
        target <- 3000
        
        # Calculate cells to remove 
        z <- sum(meta$caa_merged_total == paste0("A", region)) - target
        
        # Calculate target cells per donor
        target_per_donor <- round(target/length(unique(meta$sample_merged[meta$caa_merged_total == paste0("A", region)])))
        
        # Calculate deviation from target per donor
        caa_summary <- meta[meta$caa_merged_total == paste0("A", region),] %>% group_by(sample_merged) %>% summarize(count = n())
        caa_summary$deviation <- caa_summary$count - target_per_donor
        
        # Arrange by decreasing deviation
        caa_summary <- caa_summary %>% arrange(desc(deviation))
        
        # Filter for donors with positive deviation
        positive_summary <- caa_summary[caa_summary$deviation > 0,]
        
        # If z is large enough to reach target per donor, downsample all to target per donor
        if (z >= sum(positive_summary$deviation)) {
          positive_summary$n_sample <- target_per_donor
        } 
        
        # Else if z is large enough to reach min positive deviation, downsample to min positive deviation, then downsample equally
        else if (z >= sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation))) {
          positive_summary$n_sample <- positive_summary$count[positive_summary$deviation == min(positive_summary$deviation)]
          remaining_z <- z - sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation)) 
          if (remaining_z > 0) {
            positive_summary$n_sample <- positive_summary$n_sample - round(remaining_z/nrow(positive_summary))
          }
        }
        
        # Else if z is at least the difference between top 2 positive deviations, equalize the top 2, then downsample top 2 equally 
        else if (z >= positive_summary$deviation[1] - positive_summary$deviation[2]) {
          positive_summary$n_sample <- positive_summary$count 
          positive_summary$n_sample[1] <- positive_summary$n_sample[2] 
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
        
        # Downsample the positive donors, keep all cells for the others
        cells_keep <- c()
        for (sample in unique(meta$sample_merged[meta$caa_merged_total == paste0("A", region)])) {
          if (sample %in% positive_summary$sample_merged) {
            set.seed(100)
            cells <- sample(rownames(meta)[meta$sample_merged == sample], positive_summary$n_sample[positive_summary$sample_merged == sample], replace = FALSE)
          } else {
            cells <- rownames(meta)[meta$sample_merged == sample]
          }
          cells_keep <- c(cells_keep, cells)
        }
        
        # Keep all LCMB + downsampled CAA 
        cells_keep <- rownames(meta)[rownames(meta) %in% cells_keep | meta$sample_merged == paste0("NMA22.205.B", region)]
      }
      
      # Downsample LCMB if fold difference <= 3 but LCMB has > 3000 cells
      else if (sum(meta$sample_merged == paste0("NMA22.205.B", region)) > 3000) {
        target <- 3000
        set.seed(100)
        cells_keep <- sample(rownames(meta)[meta$sample_merged == paste0("NMA22.205.B", region)], target, replace = FALSE)
        cells_keep <- rownames(meta)[rownames(meta) %in% cells_keep | meta$caa_merged_total == paste0("A", region)] 
        
      } else { 
        cells_keep <- rownames(meta)[meta$caa_merged_total == paste0("A", region) | meta$sample_merged == paste0("NMA22.205.B", region)] 
      }
    }
    total_cells_keep <- c(total_cells_keep, cells_keep)
  }
  
  # Filter for final downsampled cells
  meta <- meta[total_cells_keep,]
  
  # Update full downsampling table
  full_meta <- rbind(full_meta, meta)
}

# Generate downsampling summary table
df <- data.frame(table(full_meta$cell_type_de, full_meta$sample_merged))
df <- df %>% pivot_wider(names_from = "Var2", values_from = "Freq") %>% data.frame()
rownames(df) <- df$Var1
df$Var1 <- NULL

df <- df[,c(paste0(c("A11.170.", "A14.193.", "NMA22.300."), "A1"), "NMA22.205.B1",
              paste0(c("A11.170.", "A14.193.", "NMA22.300."), "A3"), "NMA22.205.B3",
              paste0(c("A11.170.", "A14.193.", "NMA22.300."), "A4"), "NMA22.205.B4",
              paste0(c("A11.170.", "A14.193.", "NMA22.300."), "A9"), "NMA22.205.B9")]
write.csv(df, paste0(output_folder, "downsampled_final.csv"))

# Subset for downsampled cells 
s <- subset(s, cells = rownames(full_meta))
gc()

# Define DE group variable, make group and donor factors
s$group_de <- s$sample_merged
s$group_de[grep("\\.A1", s$group_de)] <- "A1"
s$group_de[grep("\\.A3", s$group_de)] <- "A3"
s$group_de[grep("\\.A4", s$group_de)] <- "A4"
s$group_de[grep("\\.A9", s$group_de)] <- "A9"
s$group_de[grep("\\.B1", s$group_de)] <- "B1"
s$group_de[grep("\\.B3", s$group_de)] <- "B3"
s$group_de[grep("\\.B4", s$group_de)] <- "B4"
s$group_de[grep("\\.B9", s$group_de)] <- "B9"
s$group_de <- factor(s$group_de)
s$sample_merged <- factor(s$sample_merged)

# Set idents
Idents(s) <- "group_de"

# Set default assay 
DefaultAssay(object = s) <- "SCT"

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Define comparisons to run 
comparisons <- c("B1_vs_A1", "B3_vs_A3", "B4_vs_A4", "B9_vs_A9")

# Run MAST for microglia for each comparison
for (type in c("Mg")) {

  # Subset for cell type
  s_type <- subset(s, cell_type_de == type)
  
  # Recorrect SCT data
  s_type <- PrepSCTFindMarkers(object = s_type)
  
  # Calculate CDR
  cdr <- colMeans(GetAssayData(s_type, assay = "SCT", layer = "data") > 0)
  s_type$cdr <- cdr
  
  # Standardize CDR
  s_type$cdr_centered <- scale(s_type$cdr)
  
  # Extract recorrected SCT expression data
  expressionmat_full <- GetAssayData(s_type, assay = "SCT", layer = 'data')
  
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
    LFC <- FoldChange(object = s_type, ident.1 = ident.1, ident.2 = ident.2, fc.name = "avg_log2FC", base = 2,
                      assay = "SCT", slot = "data") %>% data.frame()
    
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
    cdat <- s_type@meta.data
    cdat$wellKey <- row.names(cdat)
    fdat <- data.frame(primerid = row.names(expressionmat), row.names = row.names(expressionmat))
    
    # Create SingleCellAssay object
    sca <- FromMatrix(expressionmat, cdat, fdat, check_sanity = FALSE)
    
    # Set reference level for DE group
    group <- factor(colData(sca)$group_de)
    group <- relevel(group, ident.1) 
    colData(sca)$group_de <- group
    
    tryCatch(
      {
        # Fit model and run LRT
        zlm_group <- zlm(~ group_de + cdr_centered + (1 | sample_merged), sca, ebayes = FALSE, method = "glmer")
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
        write.csv(results, paste0(output_folder, comparison, "/results/", type, ".csv"))
      },
      error = function(e) {
        print(paste0("Error in: ", type, " ", comparison))
      }
    )
  }
}







