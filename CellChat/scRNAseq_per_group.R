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
# Summary: Infer scRNAseq cell-cell communication within AN1792, lecanemab, and nAD using CellChat
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("CellChat")
})

# Define output folder
output_folder <- "/projects/b1042/Gate_Lab/anne/AN1792_rebuttal/cellchat/scFRP_final/objects/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Load integrated scRNAseq Seurat object
s <- readRDS("/path/to/integrated/scRNAseq/object.rds")

# Define group variable
s$group <- NA
s$group[grep("^102", s$sample_merged)] <- "AN1792"
s$group[grep("\\.B", s$sample_merged)] <- "LCMB"
s$group[grep("\\.A", s$sample_merged)] <- "CAA"
print(unique(s$group))

# Define CellChat groups (broad cell types + microglia clusters) 
s$cellchat_group <- s$merged_celltype_final
s$cellchat_group[s$cellchat_group == "Mg"] <- s$fine_celltype_final[s$cellchat_group == "Mg"]

# Format cell type names
s$cellchat_group[s$cellchat_group == "Infl. Endo"] <- "Infl_Endo"
s$cellchat_group <- stringr::str_replace_all(s$cellchat_group, " ", "_")
s$cellchat_group <- stringr::str_replace_all(s$cellchat_group, "-", "_")
s$cellchat_group <- stringr::str_replace_all(s$cellchat_group, "/", "_")
print(unique(s$cellchat_group))

# Use parallel processing
future::plan("multisession", workers = 4)

# Set maximum global variable size to 1GB
options(future.globals.maxSize = 1000 * 1024^2)

# Run CellChat for each group
for (id in unique(s$group)) {
  
  # Subset for group
  s_temp <- subset(s, group == id)
  
  # Recorrect SCT data
  s_temp <- PrepSCTFindMarkers(s_temp)
  
  # Extract recorrected SCT expression and meta data
  input <- GetAssayData(s_temp, layer = "data", assay = "SCT")
  meta <- data.frame(labels = as.character(s_temp$cellchat_group), samples = factor(s_temp$sample_merged), row.names = rownames(s_temp@meta.data))
  
  # Initialize CellChat object with combined sample data
  cellchat <- createCellChat(object = input, meta = meta, group.by = "labels", datatype = "RNA")
  
  # Set CellChat database 
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB
  
  # Subset for genes in both the object and the database
  cellchat <- subsetData(cellchat)
  
  # Identify overexpressed genes and interactions
  cellchat <- identifyOverExpressedGenes(cellchat, group.by = "labels")
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Calculate communication probability
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, population.size = TRUE)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # Save object
  saveRDS(cellchat, paste0(output_folder, id, ".rds"))
}

# Run netAnalysis_computeCentrality without parallel processing 
future::plan("sequential")
for (file in list.files(output_folder)) {
  cellchat <- readRDS(paste0(output_folder, file))
  cellchat <- netAnalysis_computeCentrality(cellchat)
  saveRDS(cellchat, paste0(output_folder, file))
}


