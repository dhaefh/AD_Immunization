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
# Summary: Differential expression with MAST for plaque cluster 6
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
output_folder <- "/path/to/plaque6/mast/output/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Define filter operator
`%notin%` <- Negate(`%in%`)

# Load all-cohorts gray matter amyloid-rich integrated Seurat object and subset for cluster 6
s <- readRDS("/path/to/integrated/all/cohorts/amyloid/object.rds")
s <- subset(s, harmony_snn_res.0.4 == 6)
gc()

# Load cohort 5/7/8 data
cohort578 <- readRDS("/path/to/integrated/cohort578/object.rds")
cohort578 <- cohort578@meta.data
gc()

# Load cohort 1 data
cohort1 <- readRDS("/path/to/integrated/cohort1/object.rds")
cohort1 <- cohort1@meta.data
gc()

# Define age and sex variables for cohort 5/7/8
cohort578$age <- NA
cohort578$sex <- NA
cohort578$age[grep("^NMA22.B", cohort578$sample_id)] <- 65
cohort578$sex[grep("^NMA22.B", cohort578$sample_id)] <- "f"
cohort578$age[grep("^NMA22.A", cohort578$sample_id)] <- 82
cohort578$sex[grep("^NMA22.A", cohort578$sample_id)] <- "m"
cohort578$age[grep("^A11.170", cohort578$sample_id)] <- 66
cohort578$sex[grep("^A11.170", cohort578$sample_id)] <- "f"
cohort578$age[grep("^A14.193", cohort578$sample_id)] <- 64
cohort578$sex[grep("^A14.193", cohort578$sample_id)] <- "f"

# Define gDNA % variable for cohort 5/7/8
cohort578$gDNA_percent <- NA
cohort578$gDNA_percent[cohort578$sample_id == "NMA22.A1"] <- 1.4
cohort578$gDNA_percent[cohort578$sample_id == "NMA22.A3"] <- 1.3
cohort578$gDNA_percent[cohort578$sample_id == "NMA22.A4"] <- 2.5
cohort578$gDNA_percent[cohort578$sample_id == "NMA22.A9"] <- 1.8
cohort578$gDNA_percent[cohort578$sample_id == "NMA22.B1"] <- 1.0
cohort578$gDNA_percent[cohort578$sample_id == "NMA22.B3"] <- 0.1
cohort578$gDNA_percent[cohort578$sample_id == "NMA22.B4"] <- 0.2
cohort578$gDNA_percent[cohort578$sample_id == "NMA22.B9"] <- 0.5
cohort578$gDNA_percent[cohort578$sample_id == "A14.193.1"] <- 15.3
cohort578$gDNA_percent[cohort578$sample_id == "A14.193.3"] <- 16.3
cohort578$gDNA_percent[cohort578$sample_id == "A14.193.4"] <- 10.1
cohort578$gDNA_percent[cohort578$sample_id == "A14.193.9"] <- 1.8
cohort578$gDNA_percent[cohort578$sample_id == "A11.170.1"] <- 4.4
cohort578$gDNA_percent[cohort578$sample_id == "A11.170.3"] <- 4.9
cohort578$gDNA_percent[cohort578$sample_id == "A11.170.4"] <- 0.6
cohort578$gDNA_percent[cohort578$sample_id == "A11.170.9"] <- 0.8

# Add meta data to full Seurat object
cohort578 <- cohort578[rownames(cohort578) %in% rownames(s@meta.data),]
cohort1 <- cohort1[rownames(cohort1) %in% rownames(s@meta.data),]
print(nrow(cohort1) + nrow(cohort578) == nrow(s@meta.data))
cohort1 <- cohort1[,c("sample_barcode", "amyloid_fluo", "age", "sex", "gDNA_percent")]
colnames(cohort1)[2] <- "amyloid"
cohort578 <- cohort578[,c("sample_barcode", "cortical_amyloid", "age", "sex", "gDNA_percent")]
colnames(cohort578)[2] <- "amyloid"
meta <- rbind(cohort1, cohort578)
meta <- meta[rownames(s@meta.data),]
s$amyloid <- meta$amyloid
s$sex <- meta$sex
s$age <- meta$age
s$gDNA_percent <- meta$gDNA_percent

# Downsample CAA and LCMB spots with highest amyloid density

# Ensure no CAA donor makes up more than 50% of the CAA group
s$caa_merged <- s$sample_id
s$caa_merged[s$caa_merged %in% c("A11.170.1", "A14.193.1")] <- "A1_new"
s$caa_merged[s$caa_merged %in% c("A11.170.3", "A14.193.3")] <- "A3_new"
s$caa_merged[s$caa_merged %in% c("A11.170.4", "A14.193.4")] <- "A4_new"
s$caa_merged[s$caa_merged %in% c("A11.170.9", "A14.193.9")] <- "A9_new"

meta <- s@meta.data[s$condition %in% c("LCMB", "CAA"),]
print(unique(meta[,c("sample_id", "caa_merged")]))
old_caa <- c()
new_caa <- c()
for (region in c(1, 3, 4, 9)) {
  
  # If old CAA donor has more cells, downsample to total new controls; if new donors combined have more cells, ensure one donor does not make up more than 50% of the pool
  if (sum(meta$caa_merged == paste0("NMA22.A", region)) > sum(meta$caa_merged == paste0("A", region, "_new"))) {
    cur_meta <- meta[meta$caa_merged == paste0("NMA22.A", region),]
    cur_meta <- cur_meta %>% dplyr::arrange(desc(amyloid))
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
      cur_meta <- cur_meta %>% dplyr::arrange(desc(amyloid))
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
print(unique(meta[,c("sample_id", "caa_merged", "caa_merged_total")]))

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
        cur_meta <- cur_meta %>% dplyr::arrange(desc(amyloid))
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
    cur_meta <- cur_meta %>% dplyr::arrange(desc(amyloid))
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
        cur_meta <- cur_meta %>% dplyr::arrange(desc(amyloid))
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


# Generate table of spots per donor, ordered by region 
df <- data.frame(table(meta$sample_id))
colnames(df) <- c("sample_id", "downsampled")
rownames(df) <- df$sample_id
df <- df[c(paste0(c("A11.170.", "A14.193."), "1"), "NMA22.A1", "NMA22.B1",
           paste0(c("A11.170.", "A14.193."), "3"), "NMA22.A3", "NMA22.B3",
           paste0(c("A11.170.", "A14.193."), "4"), "NMA22.A4", "NMA22.B4",
           paste0(c("A11.170.", "A14.193."), "9"), "NMA22.A9", "NMA22.B9"),]
df$sample_id <- NULL

df2 <- data.frame(table(s$sample_id[s$condition %in% c("CAA", "LCMB")]))
colnames(df2) <- c("sample_id", "raw")
rownames(df2) <- df2$sample_id
df2 <- df2[rownames(df),]
df$raw <- df2$raw
write.csv(df[,c("raw", "downsampled")], paste0(output_folder, "lcmb_caa_downsampling.csv"))

# Manually downsample remaining donors (iAD and nAD)

# Downsample A34933.2 to make up 50% of nAD group
sample_meta <- s@meta.data[s$sample_id == "A34933.2",]
n_spots <- sum(s$condition == "nAD" & s$sample_id != "A34933.2")
sample_meta <- sample_meta %>% dplyr::arrange(desc(amyloid))
cells_keep <- rownames(sample_meta)[1:n_spots]
cells_keep <- rownames(s@meta.data)[(s$condition %in% c("iAD", "nAD") & s$sample_id != "A34933.2") | rownames(s@meta.data) %in% cells_keep]
total_cells_keep <- c(total_cells_keep, cells_keep)

# Generate table of spots per donor, ordered by condition 
temp <- s@meta.data[rownames(s@meta.data) %in% total_cells_keep & s$condition %in% c("iAD", "nAD"),]
df <- data.frame(table(temp$sample_id))
colnames(df) <- c("sample_id", "downsampled")
rownames(df) <- df$sample_id
df <- df[c(df$sample_id[-grep("^AN1792", df$sample_id)], df$sample_id[grep("^AN1792", df$sample_id)]),]
df$sample_id <- NULL

df2 <- data.frame(table(s$sample_id[s$condition %in% c("iAD", "nAD")]))
colnames(df2) <- c("sample_id", "raw")
rownames(df2) <- df2$sample_id
df2 <- df2[rownames(df),]
df$raw <- df2$raw
write.csv(df[,c("raw", "downsampled")], paste0(output_folder, "iAD_nAD_downsampling.csv"))

# Subset for downsampled spots
s <- subset(s, cells = total_cells_keep)
gc()

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5)

# Set default assay
DefaultAssay(s) <- "SCT"
for (name in names(s@assays$SCT@SCTModel.list)) {
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}

# Define region variable
s$region <- "FCX" # nAD and iAD are all FCX
s$region[s$sample_id %in% paste0(c("A14.193.", "A11.170.", "NMA22.A", "NMA22.B"), 1)] <- "FCX"
s$region[s$sample_id %in% paste0(c("A14.193.", "A11.170.", "NMA22.A", "NMA22.B"), 3)] <- "TCX"
s$region[s$sample_id %in% paste0(c("A14.193.", "A11.170.", "NMA22.A", "NMA22.B"), 4)] <- "PCX"
s$region[s$sample_id %in% paste0(c("A14.193.", "A11.170.", "NMA22.A", "NMA22.B"), 9)] <- "HIPP"
print(unique(s@meta.data[,c("sample_id", "region")]))

# Run MAST for each comparison
for (comparison in c("iAD_vs_nAD", "LCMB_vs_CAA")) {
  
  # Subset by cohort 
  if (comparison == "iAD_vs_nAD") {
    cur_s <- subset(s, condition %in% c("iAD", "nAD"))
  } else {
    cur_s <- subset(s, condition %in% c("LCMB", "CAA"))
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
  cur_s$manual_layer <- factor(cur_s$manual_layer)
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
        write.csv(results, paste0(output_folder, comparison, "/results/plaque_6.csv"))
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
        write.csv(results, paste0(output_folder, comparison, "/results/plaque_6.csv"))
      },
      error = function(e) {
        print(paste0("Error in ", comparison))
      }
    )
  }
}




  
  
  
  
  
