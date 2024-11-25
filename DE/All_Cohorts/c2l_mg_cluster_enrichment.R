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
# Summary: Differential expression with MAST for C2L enriched microglia clusters in Aß plaque niche
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
output_folder <- "/path/to/mg/cluster/mast/output/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Load all-cohorts integrated Seurat object and subset for gray matter 
s <- readRDS("/path/to/integrated/all/cohorts/object.rds")
gray_layers <- unique(s$manual_layer[grep("^gray", s$manual_layer)])
s <- subset(s, manual_layer %in% gray_layers)
gc()

# Subset for Aß-rich spots and first + second order neighbors
s <- subset(s, amyloid_neighbor_final %in% c("amyloid", "first_neighbor", "second_neighbor"))
gc()

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5)

# Create variable merging new CAA donors for downsampling
s$caa_merged <- s$sample_id
s$caa_merged[s$caa_merged %in% c("A11.170.1", "A14.193.1")] <- "A1_new"
s$caa_merged[s$caa_merged %in% c("A11.170.3", "A14.193.3")] <- "A3_new"
s$caa_merged[s$caa_merged %in% c("A11.170.4", "A14.193.4")] <- "A4_new"
s$caa_merged[s$caa_merged %in% c("A11.170.9", "A14.193.9")] <- "A9_new"
full_meta <- s@meta.data[s$condition %in% c("LCMB", "CAA"),]
print(unique(full_meta[,c("sample_id", "caa_merged")]))

# Downsample LCMB and CAA spots enriched for microglia clusters 2 and 4

# Mg-2
mg2 <- full_meta[full_meta$Mg.2_enriched == 1,]

# Ensure no CAA donor makes up more than 50% of the CAA group within each region
old_caa <- c()
new_caa <- c()
for (region in c(1, 3, 4, 9)) {
  
  # If old CAA donor has more spots, downsample to total new controls, otherwise ensure no new donors make up more than 50% of the total CAA pool
  if (sum(mg2$caa_merged == paste0("NMA22.A", region)) > sum(mg2$caa_merged == paste0("A", region, "_new"))) {
    cur_meta <- mg2[mg2$caa_merged == paste0("NMA22.A", region),]
    cur_downsample <- sum(mg2$caa_merged == paste0("A", region, "_new"))
    set.seed(100)
    cells <- sample(rownames(cur_meta), cur_downsample, replace = FALSE)
    old_caa <- c(old_caa, cells)
    new_caa <- c(new_caa, rownames(mg2)[mg2$caa_merged == paste0("A", region, "_new")])
  } else {
    
    # Calculate spots per new CAA donor
    total_pool <- sum(mg2$caa_merged %in% c(paste0("A", region, "_new"), paste0("NMA22.A", region))) 
    caa_summary <- mg2[mg2$caa_merged == paste0("A", region, "_new"),] %>% group_by(sample_id) %>% summarize(count = n()) 
    
    # Downsample if one donor makes up more than 50% (rounding up to exclude border cases) of the pool
    if (max(caa_summary$count) <= round(ceiling(total_pool/2))) {
      new_caa <- c(new_caa, rownames(mg2)[mg2$caa_merged == paste0("A", region, "_new")])
      old_caa <- c(old_caa, rownames(mg2)[mg2$caa_merged == paste0("NMA22.A", region)])
    } else {
      sample <- caa_summary$sample_id[caa_summary$count > round(ceiling(total_pool/2))]
      n_cells <- sum(mg2$caa_merged == paste0("NMA22.A", region) | (mg2$caa_merged == paste0("A", region, "_new") & mg2$sample_id != sample)) 
      cur_meta <- mg2[mg2$sample_id == sample,]
      set.seed(100)
      cells <- sample(rownames(cur_meta), n_cells, replace = FALSE)
      new_caa <- c(new_caa, rownames(mg2)[rownames(mg2) %in% cells | (mg2$caa_merged == paste0("A", region, "_new") & mg2$sample_id != sample)])
      old_caa <- c(old_caa, rownames(mg2)[mg2$caa_merged == paste0("NMA22.A", region)])
    }
  }
}

# Filter for downsampled controls  
mg2_downsampled <- mg2[rownames(mg2) %in% c(old_caa, new_caa) | str_detect(mg2$sample_id, "B"),]

# Define variable for total CAA per region
mg2_downsampled$caa_merged_total <- mg2_downsampled$caa_merged
mg2_downsampled$caa_merged_total[grep("A1", mg2_downsampled$caa_merged_total)] <- "A1"
mg2_downsampled$caa_merged_total[grep("A3", mg2_downsampled$caa_merged_total)] <- "A3"
mg2_downsampled$caa_merged_total[grep("A4", mg2_downsampled$caa_merged_total)] <- "A4"
mg2_downsampled$caa_merged_total[grep("A9", mg2_downsampled$caa_merged_total)] <- "A9"
print(unique(mg2_downsampled[,c("sample_id", "caa_merged", "caa_merged_total")]))

# Within each region, ensure fold difference < 3 between groups and each group has no more than 3000 spots
total_cells_keep <- c()
for (region in c(1, 3, 4, 9)) {
  
  # Identify larger group
  if (sum(mg2_downsampled$caa_merged_total == paste0("A", region)) > sum(mg2_downsampled$sample_id == paste0("NMA22.B", region))) {
    max_group <- "CAA"
    max <- sum(mg2_downsampled$caa_merged_total == paste0("A", region))
    min <- sum(mg2_downsampled$sample_id == paste0("NMA22.B", region))
  } else {
    max_group <- "LCMB"
    max <- sum(mg2_downsampled$sample_id == paste0("NMA22.B", region))
    min <- sum(mg2_downsampled$caa_merged_total == paste0("A", region))
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
    
    # Calculate spots to remove 
    z <- max - target
    
    # Calculate target spots per donor
    target_per_donor <- round(target/length(unique(mg2_downsampled$sample_id[mg2_downsampled$caa_merged_total == paste0("A", region)])))
    
    # Calculate deviation from target per donor
    caa_summary <- mg2_downsampled[mg2_downsampled$caa_merged_total == paste0("A", region),] %>% group_by(sample_id) %>% summarize(count = n())
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
    
    # Downsample the positive donors, keep all spots for the others
    cells_keep <- c()
    for (sample in unique(mg2_downsampled$sample_id[mg2_downsampled$caa_merged_total == paste0("A", region)])) {
      if (sample %in% positive_summary$sample_id) {
        
        cur_downsample <- positive_summary$n_sample[positive_summary$sample_id == sample]
        cur_meta <- mg2_downsampled[mg2_downsampled$sample_id == sample,]
        set.seed(100)
        cells <- sample(rownames(cur_meta), cur_downsample, replace = FALSE)
        
      } else {
        cells <- rownames(mg2_downsampled)[mg2_downsampled$sample_id == sample]
      }
      cells_keep <- c(cells_keep, cells)
    }
    
    # Keep all LCMB + downsampled CAA 
    cells_keep <- rownames(mg2_downsampled)[rownames(mg2_downsampled) %in% cells_keep | mg2_downsampled$sample_id == paste0("NMA22.B", region)]
  }
  # If fold difference > 3 and LCMB is the larger group, downsample LCMB to 3x CAA
  else if (fold > 3 & max_group == "LCMB") { 
    target <- min*3
    
    # Adjust target if greater than 3000
    if (target > 3000) {
      target <- 3000
    }
    cur_meta <- mg2_downsampled[mg2_downsampled$sample_id == paste0("NMA22.B", region),]
    set.seed(100)
    cells_keep <- sample(rownames(cur_meta), target, replace = FALSE)
    cells_keep <- rownames(mg2_downsampled)[rownames(mg2_downsampled) %in% cells_keep | mg2_downsampled$caa_merged_total == paste0("A", region)] 
  } else { 
    cells_keep <- rownames(mg2_downsampled)[mg2_downsampled$caa_merged_total == paste0("A", region) | mg2_downsampled$sample_id == paste0("NMA22.B", region)] 
  }
  total_cells_keep <- c(total_cells_keep, cells_keep)
}

# Filter for final downsampled spots
mg2_downsampled <- mg2_downsampled[total_cells_keep,]

raw <- data.frame(table(mg2$sample_id))
rownames(raw) <- raw$Var1
raw <- raw[c(paste0(c("A11.170.", "A14.193."), "1"), "NMA22.A1", "NMA22.B1",
             paste0(c("A11.170.", "A14.193."), "3"), "NMA22.A3", "NMA22.B3",
             paste0(c("A11.170.", "A14.193."), "4"), "NMA22.A4", "NMA22.B4",
             paste0(c("A11.170.", "A14.193."), "9"), "NMA22.A9", "NMA22.B9"),]

downsampled <- data.frame(table(mg2_downsampled$sample_id))
rownames(downsampled) <- downsampled$Var1
downsampled <- downsampled[c(paste0(c("A11.170.", "A14.193."), "1"), "NMA22.A1", "NMA22.B1",
                             paste0(c("A11.170.", "A14.193."), "3"), "NMA22.A3", "NMA22.B3",
                             paste0(c("A11.170.", "A14.193."), "4"), "NMA22.A4", "NMA22.B4",
                             paste0(c("A11.170.", "A14.193."), "9"), "NMA22.A9", "NMA22.B9"),]

df <- data.frame(sample_id = raw$Var1, raw = raw$Freq, downsampled = downsampled$Freq)
write.csv(df, paste0(output_folder, "mg2_cohort578_summary.csv"), row.names = FALSE)

# Mg-4
mg4 <- full_meta[full_meta$Mg.4_enriched == 1,]

# Remove samples with only one spot (causes error in PrepSCTFindMarkers) 
summary <- data.frame(table(mg4$sample_id))
samples_keep <- as.character(summary$Var1[summary$Freq >= 2])
mg4 <- mg4[mg4$sample_id %in% samples_keep,]

# Ensure no CAA donor makes up more than 50% of the CAA group within each region
old_caa <- c()
new_caa <- c()
for (region in c(1, 3, 4, 9)) {
  
  # If old CAA donor has more spots, downsample to total new controls, otherwise ensure no new donors make up more than 50% of the total CAA pool
  if (sum(mg4$caa_merged == paste0("NMA22.A", region)) > sum(mg4$caa_merged == paste0("A", region, "_new"))) {
    cur_meta <- mg4[mg4$caa_merged == paste0("NMA22.A", region),]
    cur_downsample <- sum(mg4$caa_merged == paste0("A", region, "_new"))
    set.seed(100)
    cells <- sample(rownames(cur_meta), cur_downsample, replace = FALSE)
    old_caa <- c(old_caa, cells)
    new_caa <- c(new_caa, rownames(mg4)[mg4$caa_merged == paste0("A", region, "_new")])
  } else {
    
    # Calculate spots per new CAA donor
    total_pool <- sum(mg4$caa_merged %in% c(paste0("A", region, "_new"), paste0("NMA22.A", region))) 
    caa_summary <- mg4[mg4$caa_merged == paste0("A", region, "_new"),] %>% group_by(sample_id) %>% summarize(count = n()) 
    
    # Downsample if one donor makes up more than 50% (rounding up to exclude border cases) of the pool
    if (max(caa_summary$count) <= round(ceiling(total_pool/2))) {
      new_caa <- c(new_caa, rownames(mg4)[mg4$caa_merged == paste0("A", region, "_new")])
      old_caa <- c(old_caa, rownames(mg4)[mg4$caa_merged == paste0("NMA22.A", region)])
    } else {
      sample <- caa_summary$sample_id[caa_summary$count > round(ceiling(total_pool/2))]
      n_cells <- sum(mg4$caa_merged == paste0("NMA22.A", region) | (mg4$caa_merged == paste0("A", region, "_new") & mg4$sample_id != sample)) 
      cur_meta <- mg4[mg4$sample_id == sample,]
      set.seed(100)
      cells <- sample(rownames(cur_meta), n_cells, replace = FALSE)
      new_caa <- c(new_caa, rownames(mg4)[rownames(mg4) %in% cells | (mg4$caa_merged == paste0("A", region, "_new") & mg4$sample_id != sample)])
      old_caa <- c(old_caa, rownames(mg4)[mg4$caa_merged == paste0("NMA22.A", region)])
    }
  }
}

# Filter for downsampled controls 
mg4_downsampled <- mg4[rownames(mg4) %in% c(old_caa, new_caa) | str_detect(mg4$sample_id, "B"),]

# Define variable for total CAA per region
mg4_downsampled$caa_merged_total <- mg4_downsampled$caa_merged
mg4_downsampled$caa_merged_total[grep("A1", mg4_downsampled$caa_merged_total)] <- "A1"
mg4_downsampled$caa_merged_total[grep("A3", mg4_downsampled$caa_merged_total)] <- "A3"
mg4_downsampled$caa_merged_total[grep("A4", mg4_downsampled$caa_merged_total)] <- "A4"
mg4_downsampled$caa_merged_total[grep("A9", mg4_downsampled$caa_merged_total)] <- "A9"
print(unique(mg4_downsampled[,c("sample_id", "caa_merged", "caa_merged_total")]))

# Within each region, ensure fold difference < 3 between groups and each group has no more than 3000 spots
total_cells_keep <- c()
for (region in c(1, 3, 4, 9)) {
  
  # Identify larger group
  if (sum(mg4_downsampled$caa_merged_total == paste0("A", region)) > sum(mg4_downsampled$sample_id == paste0("NMA22.B", region))) {
    max_group <- "CAA"
    max <- sum(mg4_downsampled$caa_merged_total == paste0("A", region))
    min <- sum(mg4_downsampled$sample_id == paste0("NMA22.B", region))
  } else {
    max_group <- "LCMB"
    max <- sum(mg4_downsampled$sample_id == paste0("NMA22.B", region))
    min <- sum(mg4_downsampled$caa_merged_total == paste0("A", region))
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
    
    # Calculate spots to remove 
    z <- max - target
    
    # Calculate target spots per donor
    target_per_donor <- round(target/length(unique(mg4_downsampled$sample_id[mg4_downsampled$caa_merged_total == paste0("A", region)])))
    
    # Calculate deviation from target per donor
    caa_summary <- mg4_downsampled[mg4_downsampled$caa_merged_total == paste0("A", region),] %>% group_by(sample_id) %>% summarize(count = n())
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
    
    # Downsample the positive donors, keep all spots for the others
    cells_keep <- c()
    for (sample in unique(mg4_downsampled$sample_id[mg4_downsampled$caa_merged_total == paste0("A", region)])) {
      if (sample %in% positive_summary$sample_id) {
        
        cur_downsample <- positive_summary$n_sample[positive_summary$sample_id == sample]
        cur_meta <- mg4_downsampled[mg4_downsampled$sample_id == sample,]
        set.seed(100)
        cells <- sample(rownames(cur_meta), cur_downsample, replace = FALSE)
        
      } else {
        cells <- rownames(mg4_downsampled)[mg4_downsampled$sample_id == sample]
      }
      cells_keep <- c(cells_keep, cells)
    }
    
    # Keep all LCMB + downsampled CAA 
    cells_keep <- rownames(mg4_downsampled)[rownames(mg4_downsampled) %in% cells_keep | mg4_downsampled$sample_id == paste0("NMA22.B", region)]
  }
  # If fold difference > 3 and LCMB is the larger group, downsample LCMB to 3x CAA
  else if (fold > 3 & max_group == "LCMB") { 
    target <- min*3
    
    # Adjust target if greater than 3000
    if (target > 3000) {
      target <- 3000
    }
    cur_meta <- mg4_downsampled[mg4_downsampled$sample_id == paste0("NMA22.B", region),]
    set.seed(100)
    cells_keep <- sample(rownames(cur_meta), target, replace = FALSE)
    cells_keep <- rownames(mg4_downsampled)[rownames(mg4_downsampled) %in% cells_keep | mg4_downsampled$caa_merged_total == paste0("A", region)] 
  } else { 
    cells_keep <- rownames(mg4_downsampled)[mg4_downsampled$caa_merged_total == paste0("A", region) | mg4_downsampled$sample_id == paste0("NMA22.B", region)] 
  }
  total_cells_keep <- c(total_cells_keep, cells_keep)
}

# Filter for final downsampled spots
mg4_downsampled <- mg4_downsampled[total_cells_keep,]

raw <- data.frame(table(mg4$sample_id))
rownames(raw) <- raw$Var1
raw <- raw[c(paste0(c("A11.170.", "A14.193."), "1"), "NMA22.A1", "NMA22.B1",
             paste0(c("A11.170.", "A14.193."), "3"), "NMA22.A3", "NMA22.B3",
             paste0(c("A11.170.", "A14.193."), "4"), "NMA22.A4", "NMA22.B4",
             paste0(c("A11.170.", "A14.193."), "9"), "NMA22.A9", "NMA22.B9"),] 

downsampled <- data.frame(table(mg4_downsampled$sample_id))
rownames(downsampled) <- downsampled$Var1
downsampled <- downsampled[c(paste0(c("A11.170.", "A14.193."), "1"), "NMA22.A1", "NMA22.B1",
                             paste0(c("A11.170.", "A14.193."), "3"), "NMA22.A3", "NMA22.B3",
                             paste0(c("A11.170.", "A14.193."), "4"), "NMA22.A4", "NMA22.B4",
                             paste0(c("A11.170.", "A14.193."), "9"), "NMA22.A9", "NMA22.B9"),] 

df <- data.frame(sample_id = raw$Var1, raw = raw$Freq, downsampled = downsampled$Freq)
write.csv(df, paste0(output_folder, "mg4_cohort578_summary.csv"), row.names = FALSE)

# Export iAD and nAD summary data - no need to downsample these groups
df <- s@meta.data[s$Mg.2_enriched == 1 & s$condition %in% c("iAD", "nAD"),]
df <- df %>% dplyr::group_by(condition, sample_id) %>% dplyr::summarize(n_spots = n())
df$percent_of_group <- NA
for (condition in unique(df$condition)) {
  df$percent_of_group[df$condition == condition] <- df$n_spots[df$condition == condition]/sum(df$n_spots[df$condition == condition])
}
write.csv(df, paste0(output_folder, "mg2_cohort1_summary.csv"), row.names = FALSE)

df <- s@meta.data[s$Mg.4_enriched == 1 & s$condition %in% c("iAD", "nAD"),]
df <- df %>% dplyr::group_by(condition, sample_id) %>% dplyr::summarize(n_spots = n())
df$percent_of_group <- NA
for (condition in unique(df$condition)) {
  df$percent_of_group[df$condition == condition] <- df$n_spots[df$condition == condition]/sum(df$n_spots[df$condition == condition])
}
write.csv(df, paste0(output_folder, "mg4_cohort1_summary.csv"), row.names = FALSE)

# Generate list of spots to use for DE in each cohort
cohort578_spots <- list(Mg_2 = rownames(mg2_downsampled), Mg_4 = rownames(mg4_downsampled))
cohort1_spots <- list(Mg_2 = rownames(s@meta.data)[s$condition %in% c("iAD", "nAD") & s$Mg.2_enriched == 1],
                      Mg_4 = rownames(s@meta.data)[s$condition %in% c("iAD", "nAD") & s$Mg.4_enriched == 1])

# Set default assay
DefaultAssay(s) <- "SCT"
for (name in names(s@assays$SCT@SCTModel.list)) {
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}

# Run MAST for each cluster for each comparison
for (cluster in c("Mg_2", "Mg_4")) {
  
  for (comparison in c("iAD_vs_nAD", "LCMB_vs_CAA")) {
    
    # Subset for cohort-specific enriched spots
    if (comparison == "iAD_vs_nAD") {
      cur_s <- subset(s, cells = cohort1_spots[[cluster]])
      print(unique(cur_s$condition))
    } else {
      cur_s <- subset(s, cells = cohort578_spots[[cluster]])
      print(unique(cur_s$condition))
    }
    
    # Recorrect SCT data
    skip_to_next <- FALSE
    tryCatch(
      {
        cur_s <- PrepSCTFindMarkers(cur_s)
      },
      error = function(e) {
        print(paste0("Cannot recorrect SCT data for ", comparison))
        skip_to_next <<- TRUE
      }
    )
    if (skip_to_next) {
      print(paste0("Skipping ", comparison, ": cannot recorrect SCT data"))
      next 
    }
    
    # Calculate CDR 
    cdr <- colMeans(GetAssayData(cur_s, assay = "SCT", layer = "data") > 0)
    cur_s$cdr <- cdr
    
    # Standardize continuous covariates and make categorical covariates factors
    cur_s$gDNA_centered <- scale(cur_s$gDNA_percent)
    cur_s$age_centered <- scale(cur_s$age)
    cur_s$cdr_centered <- scale(cur_s$cdr)
    cur_s$sex <- factor(cur_s$sex)
    cur_s$sample_id <- factor(cur_s$sample_id)
    cur_s$region <- factor(cur_s$region)
    
    # Extract recorrected SCT expression data
    expressionmat_full <- GetAssayData(cur_s, assay = "SCT", layer = 'data')
    
    # Run MAST with cohort-specific covariates
    if (comparison == "iAD_vs_nAD") {
      
      # Create subfolders
      dir.create(paste0(output_folder, comparison), showWarnings = FALSE, recursive = TRUE)
      dir.create(paste0(output_folder, comparison, "/results"), showWarnings = FALSE, recursive = TRUE)
      dir.create(paste0(output_folder, comparison, "/plots"), showWarnings = FALSE, recursive = TRUE)
      
      # Extract names of comparison groups
      ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
      ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
      
      # Set idents
      Idents(cur_s) <- "condition"
      
      # Calculate fold change and percent expression using SCT data
      LFC <- FoldChange(object = cur_s, ident.1 = ident.1, ident.2 = ident.2, fc.name = "avg_log2FC", base = 2,
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
      cdat <- cur_s@meta.data
      cdat$wellKey <- row.names(cdat)
      fdat <- data.frame(primerid = row.names(expressionmat), row.names = row.names(expressionmat))
      
      # Create SingleCellAssay object
      sca <- FromMatrix(expressionmat, cdat, fdat, check_sanity = FALSE)
      
      # Set reference level for condition
      cond <- factor(colData(sca)$condition)
      cond <- relevel(cond, ident.1)
      colData(sca)$condition <- cond
      
      tryCatch(
        {
          # Fit model and run LRT
          zlm <- zlm(~ condition + sex + age_centered + cdr_centered + gDNA_centered + (1 | sample_id), sca, ebayes = FALSE, method = "glmer")
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
          write.csv(results, paste0(output_folder, comparison, "/results/", cluster, ".csv"))
        },
        error = function(e) {
          print(paste0("Error in ", comparison))
        }
      )
      
    } else {
      # Create subfolders
      dir.create(paste0(output_folder, comparison), showWarnings = FALSE, recursive = TRUE)
      dir.create(paste0(output_folder, comparison, "/results"), showWarnings = FALSE, recursive = TRUE)
      dir.create(paste0(output_folder, comparison, "/plots"), showWarnings = FALSE, recursive = TRUE)
      
      # Extract names of comparison groups
      ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
      ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
      
      # Set idents
      Idents(cur_s) <- "condition"
      
      # Calculate fold change and percent expression using SCT data
      LFC <- FoldChange(object = cur_s, ident.1 = ident.1, ident.2 = ident.2, fc.name = "avg_log2FC", base = 2,
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
      cdat <- cur_s@meta.data
      cdat$wellKey <- row.names(cdat)
      fdat <- data.frame(primerid = row.names(expressionmat), row.names = row.names(expressionmat))
      
      # Create SingleCellAssay object
      sca <- FromMatrix(expressionmat, cdat, fdat, check_sanity = FALSE)
      
      # Set reference level for condition
      cond <- factor(colData(sca)$condition)
      cond <- relevel(cond, ident.1)
      colData(sca)$condition <- cond
      
      tryCatch(
        {
          # Fit model and run LRT 
          zlm <- zlm(~ condition + region + cdr_centered + gDNA_centered + (1 | sample_id), sca, ebayes = FALSE, method = "glmer")
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
          write.csv(results, paste0(output_folder, comparison, "/results/", cluster, ".csv"))
        },
        error = function(e) {
          print(paste0("Error in ", comparison))
        }
      )
    }
  }
}





