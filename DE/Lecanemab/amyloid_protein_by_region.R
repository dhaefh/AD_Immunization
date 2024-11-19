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
# Summary: Differential protein expression with FindMarkers for amyloid-rich gray matter
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
output_folder <- "/path/to/amyloid/by/region/protein/output/folder/"

# Load cohort 5/7/8 protein Seurat object and integrated cohort 5/7/8 RNA Seurat object 
s <- readRDS("/path/to/cohort578/final/protein/object.rds")
temp <- readRDS("/path/to/integrated/cohort578/object.rds")

# Add meta data to protein object
sum(row.names(temp@meta.data) != row.names(s@meta.data))
s@meta.data <- cbind(s@meta.data, temp@meta.data[,c("manual_layer", "cortical_amyloid", "vascular_amyloid")])

# Define cortical amyloid enrichment
s$amyloid_rich <- "not_rich"
s$amyloid_rich[s$cortical_amyloid > 183 & s$vascular_amyloid == 0] <- "rich"

# Subset for cortical amyloid-rich spots in gray matter 
gray_layers <- unique(s$manual_layer[grep("gray", s$manual_layer)])
s <- subset(s, amyloid_rich == "rich" & manual_layer %in% gray_layers)
gc()

# Downsample spots with highest amyloid density

# Ensure no CAA donor makes up more than 50% of the CAA group
s$caa_merged <- s$sample_id
s$caa_merged[s$caa_merged %in% c("A11.170.1", "A14.193.1")] <- "A1_new"
s$caa_merged[s$caa_merged %in% c("A11.170.3", "A14.193.3")] <- "A3_new"
s$caa_merged[s$caa_merged %in% c("A11.170.4", "A14.193.4")] <- "A4_new"
s$caa_merged[s$caa_merged %in% c("A11.170.9", "A14.193.9")] <- "A9_new"
meta <- s@meta.data
old_caa <- c()
new_caa <- c()
for (region in c(1, 3, 4, 9)) {
  
  # If old CAA donor has more cells, downsample to total new controls; if new donors combined have more cells, ensure one donor does not make up more than 50% of the pool
  if (sum(meta$caa_merged == paste0("NMA22.A", region)) > sum(meta$caa_merged == paste0("A", region, "_new"))) {
    cur_meta <- meta[meta$caa_merged == paste0("NMA22.A", region),]
    cur_meta <- cur_meta %>% dplyr::arrange(desc(cortical_amyloid))
    cur_downsample <- sum(meta$caa_merged == paste0("A", region, "_new"))
    cells <- rownames(cur_meta)[1:cur_downsample]
    old_caa <- c(old_caa, cells)
    new_caa <- c(new_caa, rownames(meta)[meta$caa_merged == paste0("A", region, "_new")])
  } else {
    total_pool <- sum(meta$caa_merged %in% c(paste0("A", region, "_new"), paste0("NMA22.A", region))) # Full CAA pool for current cell type in current region
    caa_summary <- meta[meta$caa_merged == paste0("A", region, "_new"),] %>% group_by(sample_id) %>% summarize(count = n()) # Cells per new donor
    
    # If no donors make up more than 50% of the entire CAA pool, keep all cells (take ceiling to avoid downsampling border cases)
    if (max(caa_summary$count) <= round(ceiling(total_pool/2))) {
      new_caa <- c(new_caa, rownames(meta)[meta$caa_merged == paste0("A", region, "_new")])
      old_caa <- c(old_caa, rownames(meta)[meta$caa_merged == paste0("NMA22.A", region)])
    } else {
      # There can only be one sample that makes up more than 50% of the entire CAA pool 
      sample <- caa_summary$sample_id[caa_summary$count > round(ceiling(total_pool/2))]
      n_cells <- sum(meta$caa_merged == paste0("NMA22.A", region) | (meta$caa_merged == paste0("A", region, "_new") & meta$sample_id != sample)) # Calculate remaining cells in the CAA pool
      
      cur_meta <- meta[meta$sample_id == sample,]
      cur_meta <- cur_meta %>% dplyr::arrange(desc(cortical_amyloid))
      cells <- rownames(cur_meta)[1:n_cells]
      
      new_caa <- c(new_caa, rownames(meta)[rownames(meta) %in% cells | (meta$caa_merged == paste0("A", region, "_new") & meta$sample_id != sample)])
      old_caa <- c(old_caa, rownames(meta)[meta$caa_merged == paste0("NMA22.A", region)])
    }
  }
}

# Subset for downsampled controls - keep all LCMB cells in all regions 
meta <- meta[rownames(meta) %in% c(old_caa, new_caa) | str_detect(meta$sample_id, "B"),]

# Create variable for total CAA per region
meta$caa_merged_total <- meta$caa_merged
meta$caa_merged_total[grep("A1", meta$caa_merged_total)] <- "A1"
meta$caa_merged_total[grep("A3", meta$caa_merged_total)] <- "A3"
meta$caa_merged_total[grep("A4", meta$caa_merged_total)] <- "A4"
meta$caa_merged_total[grep("A9", meta$caa_merged_total)] <- "A9"

# Ensure fold difference < 3 between groups and each group has no more than 3000 spots (within each brain region)
total_cells_keep <- c()
for (region in c(1, 3, 4, 9)) {
  
  # Identify group with max cells in region
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
  
  # If fold difference is greater than 3 and CAA has more cells, downsample to adjusted target max
  if (fold > 3 & max_group == "CAA") {
    
    # Set target number of cells for CAA (LCMB*3)
    target <- min*3
    
    # Adjust if greater than 3000
    if (target > 3000) {
      target <- 3000
    }
    
    # Calculate z value 
    z <- max - target
    
    # Calculate target cells per donor
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
      remaining_z <- z - sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation)) # Subtract cells removed in previous step from total z
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
    
    # Downsample the positive donors, keep all cells for the others
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
  # If fold difference is greater than 3 and LCMB has more cells, downsample to adjusted target max
  else if (fold > 3 & max_group == "LCMB") { 
    target <- min*3
    
    # Adjust if greater than 3000
    if (target > 3000) {
      target <- 3000
    }
    cur_meta <- meta[meta$sample_id == paste0("NMA22.B", region),]
    cur_meta <- cur_meta %>% dplyr::arrange(desc(cortical_amyloid))
    cells_keep <- rownames(cur_meta)[1:target]
    cells_keep <- rownames(meta)[rownames(meta) %in% cells_keep | meta$caa_merged_total == paste0("A", region)] # Keep all CAA cells + downsampled LCMB
  } 
  
  # Cases where fold <= 3 but larger group has more than 3000 - for this cohort, only possible for CAA
  else if (sum(meta$caa_merged_total == paste0("A", region)) > 3000) {
    target <- 3000
    
    # Calculate z value for CAA 
    z <- sum(meta$caa_merged_total == paste0("A", region)) - target
    
    # Calculate target cells per donor
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
      remaining_z <- z - sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation)) # Subtract cells removed in previous step from total z
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
    
    # Downsample the positive donors, keep all cells for the others
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
  } else { # Otherwise keep all cells
    cells_keep <- rownames(meta)[meta$caa_merged_total == paste0("A", region) | meta$sample_id == paste0("NMA22.B", region)] # Keep all cells
  }
  total_cells_keep <- c(total_cells_keep, cells_keep)
}
meta <- meta[total_cells_keep,]


# Generate downsampling summary table
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
print(unique(s@meta.data[,c("group_de", "sample_id")]))

# Make categorical covariates factors
s$group_de <- factor(s$group_de)
s$sample_id <- factor(s$sample_id)
s$manual_layer <- factor(s$manual_layer)

# Set idents
Idents(s) <- "group_de"

# Combine counts across samples
s <- JoinLayers(s)

# Calculate CDR using isotype-normalized counts
cdr <- colMeans(GetAssayData(s, assay = "Protein", layer = "counts") > 0)
s@meta.data$cdr <- cdr

# Standardize CDR
s@meta.data$cdr_centered <- scale(s$cdr)

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Define comparisons to run
comparisons <- c("B1_vs_A1", "B3_vs_A3", "B4_vs_A4", "B9_vs_A9")

# Run FindMarkers with negative binomial test for each comparison 
for (comparison in comparisons) {
  
  # Create output folders
  dir.create(paste0(output_folder, comparison), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(output_folder, comparison, "/results"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(output_folder, comparison, "/plots"), showWarnings = FALSE, recursive = TRUE)
  
  # Extract names of comparison groups
  ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
  ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
  
  # Exclude manual layer as a covariate for hippocampus since there is only gray-hippocampus - causes error in FindMarkers
  if (comparison == "B9_vs_A9") {
    results <- FindMarkers(s, assay = "Protein", ident.1 = ident.1, ident.2 = ident.2, test.use = "negbinom", min.pct = 0.01, logfc.threshold = -Inf, 
                           latent.vars =  c("cdr_centered"), fc.slot = "counts", slot = "counts")
  } else {
    results <- FindMarkers(s, assay = "Protein", ident.1 = ident.1, ident.2 = ident.2, test.use = "negbinom", min.pct = 0.01, logfc.threshold = -Inf, 
                           latent.vars =  c("cdr_centered", "manual_layer"), fc.slot = "counts", slot = "counts")
  }
  
  
  # Define DEPs and save results
  results$BH <- p.adjust(results$p_val, method = "BH")
  results$antibody <- row.names(results)
  results$DE <- "Not DE"
  results$DE[results$avg_log2FC > fc_thresh & results$BH < p_thresh] <- "Upregulated"
  results$DE[results$avg_log2FC < -fc_thresh & results$BH < p_thresh] <- "Downregulated"
  results$DE_antibody <- NA
  results$DE_antibody[results$DE != "Not DE"] <- results$antibody[results$DE != "Not DE"]
  write.csv(results, paste0(output_folder, comparison, "/results/amyloid_rich.csv"))
}

