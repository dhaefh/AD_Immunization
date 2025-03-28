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

# Subset for AN1792 
cohort1 <- readRDS("/path/to/integrated/AN1792/object.rds")
cohort1 <- cohort1@meta.data
gc()
cohort1 <- cohort1[rownames(cohort1) %in% rownames(s@meta.data),]
s <- subset(s, cells = rownames(cohort1))
gc()

# Subset for Aß-rich spots and first + second order neighbors
s <- subset(s, amyloid_neighbor_final %in% c("amyloid", "first_neighbor", "second_neighbor"))
gc()

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Downsample spots enriched for microglia clusters 2 and 4 such that no sample makes up more than 50% of iAD-Ext and DE groups have less than 3 fold difference

# Mg-2 
df <- s@meta.data[s$Mg.2_enriched == 1 & s$condition %in% c("iAD", "nAD"),]

# Downsample nAD

# Calculate spots to remove
target <- 3*sum(df$condition_residual_amyloid == "ext")
z <- sum(df$condition == "nAD") - target

# Calculate target spots per donor
target_per_donor <- round(target/length(unique(df$sample_id[df$condition == "nAD"])))

# Calculate deviation from target per donor
summary <- df[df$condition == "nAD",] %>% dplyr::group_by(sample_id) %>% dplyr::summarize(count = n())
summary$deviation <- summary$count - target_per_donor

# Arrange by decreasing deviation
summary <- summary %>% dplyr::arrange(desc(deviation))

# Filter for donors with positive deviation
positive_summary <- summary[summary$deviation > 0,]

# Downsample to the minimum positive deviation, then downsample equally
positive_summary$n_sample <- positive_summary$count[positive_summary$deviation == min(positive_summary$deviation)]
remaining_z <- z - sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation)) 
if (remaining_z > 0) {
  positive_summary$n_sample <- positive_summary$n_sample - round(remaining_z/nrow(positive_summary))
}

# Downsample positive donors
cells_keep <- c()
for (sample in unique(df$sample_id[df$condition == "nAD"])) {
  if (sample %in% positive_summary$sample_id) {
    cur_downsample <- positive_summary$n_sample[positive_summary$sample_id == sample]
    cur_meta <- df[df$sample_id == sample,]
    set.seed(100)
    cells <- sample(rownames(cur_meta), cur_downsample, replace = FALSE)
  } else {
    cells <- rownames(df)[df$sample_id == sample]
  }
  cells_keep <- c(cells_keep, cells)
}

# Subset for downsampled nAD
df <- df[df$condition != "nAD" | rownames(df) %in% cells_keep,]

# Define final downsampled Mg-4 spots
mg2_keep <- rownames(df)

# Generate summary table
df <- df %>% dplyr::group_by(condition_residual_amyloid, sample_id) %>% dplyr::summarize(n_spots = n())
df$percent_of_group <- NA
for (condition in unique(df$condition_residual_amyloid)) {
  df$percent_of_group[df$condition_residual_amyloid == condition] <- df$n_spots[df$condition_residual_amyloid == condition]/sum(df$n_spots[df$condition_residual_amyloid == condition])
}
write.csv(df, paste0(output_folder, "mg2_cohort1_by_clearance_summary.csv"), row.names = FALSE)

# Mg-4
df <- s@meta.data[s$Mg.4_enriched == 1 & s$condition %in% c("iAD", "nAD"),]

# Downsample 102.6
set.seed(100)
ext_keep <- sample(rownames(df)[df$sample_id == "AN1792.102.6"], sum(df$sample_id != "AN1792.102.6" & df$condition_residual_amyloid == "ext"), replace = FALSE)
df <- df[df$sample_id != "AN1792.102.6" | rownames(df) %in% ext_keep,]

# Downsample nAD

# Calculate spots to remove
target <- 3*sum(df$condition_residual_amyloid == "ext")
z <- sum(df$condition == "nAD") - target

# Calculate target spots per donor
target_per_donor <- round(target/length(unique(df$sample_id[df$condition == "nAD"])))

# Calculate deviation from target per donor
summary <- df[df$condition == "nAD",] %>% dplyr::group_by(sample_id) %>% dplyr::summarize(count = n())
summary$deviation <- summary$count - target_per_donor

# Arrange by decreasing deviation
summary <- summary %>% dplyr::arrange(desc(deviation))

# Filter for donors with positive deviation
positive_summary <- summary[summary$deviation > 0,]

# Downsample to the minimum positive deviation, then downsample equally
positive_summary$n_sample <- positive_summary$count[positive_summary$deviation == min(positive_summary$deviation)]
remaining_z <- z - sum(positive_summary$deviation[positive_summary$deviation > min(positive_summary$deviation)] - min(positive_summary$deviation)) 
if (remaining_z > 0) {
  positive_summary$n_sample <- positive_summary$n_sample - round(remaining_z/nrow(positive_summary))
}

# Downsample positive donors
cells_keep <- c()
for (sample in unique(df$sample_id[df$condition == "nAD"])) {
  if (sample %in% positive_summary$sample_id) {
    cur_downsample <- positive_summary$n_sample[positive_summary$sample_id == sample]
    cur_meta <- df[df$sample_id == sample,]
    set.seed(100)
    cells <- sample(rownames(cur_meta), cur_downsample, replace = FALSE)
  } else {
    cells <- rownames(df)[df$sample_id == sample]
  }
  cells_keep <- c(cells_keep, cells)
}

# Subset for downsampled nAD
df <- df[df$condition != "nAD" | rownames(df) %in% cells_keep,]

# Downsample iAD-Lim

# Calculate spots to remove
target <- 3*sum(df$condition == "nAD")
z <- sum(df$condition_residual_amyloid == "lim") - target

# Calculate target cells per donor
target_per_donor <- round(target/length(unique(df$sample_id[df$condition_residual_amyloid == "lim"])))

# Calculate deviation from target per donor
summary <- df[df$condition_residual_amyloid == "lim",] %>% dplyr::group_by(sample_id) %>% dplyr::summarize(count = n())
summary$deviation <- summary$count - target_per_donor

# Arrange by decreasing deviation
summary <- summary %>% dplyr::arrange(desc(deviation))

# Filter for donors with positive deviation
positive_summary <- summary[summary$deviation > 0,]

# Downsample largest donor
positive_summary$n_sample <- positive_summary$count
positive_summary$n_sample[1] <- positive_summary$n_sample[1] - z
cells_keep <- c()
for (sample in unique(df$sample_id[df$condition_residual_amyloid == "lim"])) {
  if (sample %in% positive_summary$sample_id) {
    cur_downsample <- positive_summary$n_sample[positive_summary$sample_id == sample]
    cur_meta <- df[df$sample_id == sample,]
    set.seed(100)
    cells <- sample(rownames(cur_meta), cur_downsample, replace = FALSE)
  } else {
    cells <- rownames(df)[df$sample_id == sample]
  }
  cells_keep <- c(cells_keep, cells)
}

# Subset for downsampled iAD-Lim
df <- df[df$condition_residual_amyloid != "lim" | rownames(df) %in% cells_keep,]

# Define final downsampled Mg-4 spots
mg4_keep <- rownames(df)

# Generate summary table
df <- df %>% dplyr::group_by(condition_residual_amyloid, sample_id) %>% dplyr::summarize(n_spots = n())
df$percent_of_group <- NA
for (condition in unique(df$condition_residual_amyloid)) {
  df$percent_of_group[df$condition_residual_amyloid == condition] <- df$n_spots[df$condition_residual_amyloid == condition]/sum(df$n_spots[df$condition_residual_amyloid == condition])
}
write.csv(df, paste0(output_folder, "mg4_cohort1_by_clearance_summary.csv"), row.names = FALSE)

# Define list of DE spots
de_spots <- list(Mg_2 = mg2_keep, Mg_4 = mg4_keep)

# Prep for recorrection
DefaultAssay(s) <- "SCT"
for (name in names(s@assays$SCT@SCTModel.list)) {
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}

# Run MAST for each cluster for each comparison
for (cluster in c("Mg_2", "Mg_4")) {
  
  for (comparison in c("lim_vs_nAD", "ext_vs_nAD")) {
    
    # Subset for downsampled enriched spots
    cur_s <- subset(s, cells = de_spots[[cluster]])
    print(unique(cur_s$condition))
    
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
    
    # Extract recorrected SCT expression data
    expressionmat_full <- GetAssayData(cur_s, assay = "SCT", layer = 'data')
    
    # Create subfolders
    dir.create(paste0(output_folder, comparison), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(output_folder, comparison, "/results"), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(output_folder, comparison, "/plots"), showWarnings = FALSE, recursive = TRUE)
    
    # Extract names of comparison groups
    ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
    ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
    
    # Set idents
    Idents(cur_s) <- "condition_residual_amyloid"
    
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
    cond <- factor(colData(sca)$condition_residual_amyloid)
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
  }
}





