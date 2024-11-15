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
# Summary: Define cell type enrichment using Cell2Location predicted abundance
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages ({
  library('Seurat')
  library('glue')
  library('dplyr')
  library('stringr')
  library("Matrix")
  library("data.table")
})

# Define filter operator 
`%notin%` <- Negate(`%in%`)

# Load integrated cohort 1 Seurat and extract meta data
s <- readRDS("/path/to/integrated/cohort1/object.rds")
layers <- s@meta.data
s <- NULL
gc()

# Load C2L q05 cell type abundance predictions
meta <- read.csv("/path/to/cell2location/mapping/q05_meta.csv", row.names = 1)

# Add manual layer to meta data
print(sum(rownames(layers) != rownames(meta)))
meta$manual_layer <- layers$manual_annotation

# Define fibroblast, pericyte, peripheral immune cell, smooth muscle cell and endothelial cell enrichment in top 1% of gray/white matter or 5% of meninges
for (type in c("FB", "Pericyte", "Peripheral.Immune.Cell", "SMC", "endothelial")) {
  meta[[paste0(type, "_enriched")]] <- 0
  for (sample in unique(meta$sample_id)) {
    
    # Filter for current sample 
    cur_meta <- meta[meta$sample_id == sample,]
    
    # Gray/white
    region_meta <- cur_meta[cur_meta$manual_layer %notin% c("meninges", "gray-l1"),]
    
    # Calculate quantile for q05 predictions for current cell type
    quantile <- quantile(region_meta[[type]], 0.99)
    
    # Define enriched spots for current cell type in current sample 
    meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$sample_id == sample & meta$manual_layer %notin% c("meninges", "gray-l1")] <- 1
    
    # Gray/white
    region_meta <- cur_meta[cur_meta$manual_layer == "meninges",]
    
    # Calculate quantile for q05 predictions for current cell type
    quantile <- quantile(region_meta[[type]], 0.95)
    
    # Define enriched spots for current cell type in current sample 
    meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$sample_id == sample & meta$manual_layer == "meninges"] <- 1
  }
}

# Define interneuron enrichment in top 5% of gray matter
type <- "Interneuron"
meta[[paste0(type, "_enriched")]] <- 0
for (sample in unique(meta$sample_id)) {
  
  # Filter for current sample 
  cur_meta <- meta[meta$sample_id == sample,]
  
  # Filter for ROI 
  region_meta <- cur_meta[cur_meta$manual_layer %notin% c("white", "meninges", "gray-l1"),]
  
  # Calculate quantile for q05 predictions for current cell type
  quantile <- quantile(region_meta[[type]], 0.95)
  
  # Define enriched spots for current cell type in current sample 
  meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$manual_layer %notin% c("white", "meninges", "gray-l1") & meta$sample_id == sample] <- 1
}

# Define layer specific neuron enrichment in relevant cortical layers, defining l2/3 neurons in top 10%, l4 EN in top 50%, l4/5 EN in top 15%
for (type in c("L2.3.EN", "L4.EN", "L4.5.EN")) {
  meta[[paste0(type, "_enriched")]] <- 0
  for (sample in unique(meta$sample_id)) {
    
    # Filter for current sample 
    cur_meta <- meta[meta$sample_id == sample,]
    
    # Identify ROI and quantile
    if (type == "L2.3.EN") {
      layers <- c("gray-l2", "gray-l3")
      cutoff <- 0.9
    } else if (type == "L4.EN") {
      layers <- c("gray-l4")
      cutoff <- 0.5
    } else {
      layers <- c("gray-l4", "gray-l56")
      cutoff <- 0.85
    }
    
    # Filter for ROI 
    region_meta <- cur_meta[cur_meta$manual_layer %in% layers,]
    
    # Calculate quantile for q05 predictions for current cell type
    quantile <- quantile(region_meta[[type]], cutoff)
    
    # Define enriched spots for current cell type in current sample 
    meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$manual_layer %in% layers & meta$sample_id == sample] <- 1
  }
}

# Define l5 EN enrichment in top 5% of relevant cortical layers
type <- "L5.EN"
meta[[paste0(type, "_enriched")]] <- 0
for (sample in unique(meta$sample_id)) {
  
  # Filter for current sample 
  cur_meta <- meta[meta$sample_id == sample,]
  
  # Filter for ROI 
  region_meta <- cur_meta[cur_meta$manual_layer == "gray-l56",]
  
  # Calculate quantile for q05 predictions for current cell type
  quantile <- quantile(region_meta[[type]], 0.95)
  
  # Define enriched spots for current cell type in current sample 
  meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$manual_layer == "gray-l56" & meta$sample_id == sample] <- 1
}

# Define l5/6 EN and l5/6 CCa EN enrichment in top 5% of relevant cortical layers 
for (type in c("L5.6.EN", "L5.6.CCa.EN")) {
  meta[[paste0(type, "_enriched")]] <- 0
  for (sample in unique(meta$sample_id)) {
    
    # Filter for current sample 
    cur_meta <- meta[meta$sample_id == sample,]
    
    # Filter for ROI 
    region_meta <- cur_meta[cur_meta$manual_layer == "gray-l56",]
    
    # Calculate quantile for q05 predictions for current cell type
    quantile <- quantile(region_meta[[type]], 0.95)
    
    # Define enriched spots for current cell type in current sample 
    meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$manual_layer == "gray-l56" & meta$sample_id == sample] <- 1
  }
}

# Define l5/6 CCb EN enrichment in top 5% of relevant cortical layers 
type <- "L5.6.CCb.EN"
meta[[paste0(type, "_enriched")]] <- 0
for (sample in unique(meta$sample_id)) {
  
  # Filter for current sample 
  cur_meta <- meta[meta$sample_id == sample,]
  
  # Filter for ROI 
  region_meta <- cur_meta[cur_meta$manual_layer == "gray-l56",]
  
  # Calculate quantile for q05 predictions for current cell type
  quantile <- quantile(region_meta[[type]], 0.95)
  
  # Define enriched spots for current cell type in current sample 
  meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$manual_layer == "gray-l56" & meta$sample_id == sample] <- 1
}

# Define MYO16 EN enrichment in top 1% of gray matter
type <- "MYO16.EN"
meta[[paste0(type, "_enriched")]] <- 0
for (sample in unique(meta$sample_id)) {
  
  # Filter for current sample 
  cur_meta <- meta[meta$sample_id == sample,]
  
  # Filter for ROI 
  region_meta <- cur_meta[cur_meta$manual_layer %notin% c("white", "meninges", "gray-l1"),]
  
  # Calculate quantile for q05 predictions for current cell type
  quantile <- quantile(region_meta[[type]], 0.99)
  
  # Define enriched spots for current cell type in current sample 
  meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$manual_layer %notin% c("white", "meninges", "gray-l1") & meta$sample_id == sample] <- 1
}

# Define oligodendrocyte precursor cell enrichment in top 5% of gray/white matter
type <- "OPC"
meta[[paste0(type, "_enriched")]] <- 0
for (sample in unique(meta$sample_id)) {
  
  # Filter for current sample 
  cur_meta <- meta[meta$sample_id == sample,]
  
  # Filter for ROI 
  region_meta <- cur_meta[cur_meta$manual_layer %notin% c("meninges", "gray-l1"),]
  
  # Calculate quantile for q05 predictions for current cell type
  quantile <- quantile(region_meta[[type]], 0.95)
  
  # Define enriched spots for current cell type in current sample 
  meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$manual_layer %notin% c("meninges", "gray-l1") & meta$sample_id == sample] <- 1
}

# Define oligodendrocyte enrichment in top 30% of white matter
type <- "oligodendrocyte"
meta[[paste0(type, "_enriched")]] <- 0
for (sample in unique(meta$sample_id)) {
  
  # Filter for current sample 
  cur_meta <- meta[meta$sample_id == sample,]
  
  # Filter for ROI 
  region_meta <- cur_meta[cur_meta$manual_layer == "white",]
  
  # Calculate quantile for q05 predictions for current cell type
  quantile <- quantile(region_meta[[type]], 0.7)
  
  # Define enriched spots for current cell type in current sample 
  meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$manual_layer == "white" & meta$sample_id == sample] <- 1
}

# Define astrocyte enrichment in top 5% of gray matter or 5% of white matter
type <- "Astrocyte"
meta[[paste0(type, "_enriched")]] <- 0
for (sample in unique(meta$sample_id)) {
  
  # Filter for current sample 
  cur_meta <- meta[meta$sample_id == sample,]
  
  # Filter for white matter  
  region_meta <- cur_meta[cur_meta$manual_layer == "white",]
  
  # Calculate quantile for q05 predictions for current cell type
  quantile <- quantile(region_meta[[type]], 0.95)
  
  # Define enriched spots in white matter 
  meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$manual_layer == "white" & meta$sample_id == sample] <- 1
  
  # Filter for gray matter  
  region_meta <- cur_meta[cur_meta$manual_layer %notin% c("white", "meninges", "gray-l1"),]
  
  # Calculate quantile for q05 predictions for current cell type
  quantile <- quantile(region_meta[[type]], 0.95)
  
  # Define enriched spots in white matter 
  meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$manual_layer %notin% c("white", "meninges", "gray-l1") & meta$sample_id == sample] <- 1
}

# Define microglia enrichment in top 5% of gray matter or 5% of white matter
type <- "Microglia"
meta[[paste0(type, "_enriched")]] <- 0
for (sample in unique(meta$sample_id)) {
  
  # Filter for current sample 
  cur_meta <- meta[meta$sample_id == sample,]
  
  # Filter for white matter  
  region_meta <- cur_meta[cur_meta$manual_layer == "white",]
  
  # Calculate quantile for q05 predictions for current cell type
  quantile <- quantile(region_meta[[type]], 0.95)
  
  # Define enriched spots in white matter 
  meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$manual_layer == "white" & meta$sample_id == sample] <- 1
  
  # Filter for gray matter  
  region_meta <- cur_meta[cur_meta$manual_layer %notin% c("white", "meninges", "gray-l1"),]
  
  # Calculate quantile for q05 predictions for current cell type
  quantile <- quantile(region_meta[[type]], 0.95)
  
  # Define enriched spots in white matter 
  meta[[paste0(type, "_enriched")]][meta[[type]] > quantile & meta$manual_layer %notin% c("white", "meninges", "gray-l1") & meta$sample_id == sample] <- 1
}

# Save meta data with cell type enrichment 
write.csv(meta, "/path/to/cell2location/mapping/q05_meta_with_enrichment.csv")


