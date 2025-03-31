# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                         AN1792 Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Written by: Anne Forsyth 
# Summary: Generate UpSet plot using DESeq2 or MAST output
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ggrepel")
  library("UpSetR")
  library("ComplexUpset")
  library("randomcoloR")
  library("stringr")
})


# Define path to folder containing csv files of all results to use in upset plot
# Ensure that the only files in this folder are the results, and the names correspond to desired set names in plot
result_dir <- "/path/to/all/results/folder/"

# Set these variables to column names corresponding to LFC, adjusted p-value and gene name (so they can be renamed in the plot function)
lfc <- "<name of LFC column>"
padj <- "<name of padj column>"
gene <- "<name of gene name column>"

# Define DE thresholds
fc_thresh <- log2(1.5) 
p_thresh <- 0.05

# Define function to generate list of DEGs from results
load_degs <- function(result_dir, pct_in_title = FALSE, n_in_title = FALSE) {
  
  # Initialize list for DEG sets 
  deg_list <- list()
  
  # Load all results and add DEGs to list
  for (file in list.files(result_dir)) {
    
    # Load results
    results <- read.csv(paste0(result_dir, file), row.names = 1)
    
    # Rename variables 
    results$lfc <- results[[lfc]]
    results$padj <- results[[padj]]
    results$gene <- results[[gene]]
    
    # Define DEGs
    results$DE <- "Not DE"
    results$DE[results$padj < p_thresh & abs(results$lfc) > fc_thresh] <- "DE"
    
    # Add to upset list
    set_name <- str_replace_all(file, ".csv", "")
    degs <- results$gene[results$DE == "DE"]
    
    # Customize set name
    if (pct_in_title & n_in_title) {
      stop("Specify only one metric to include in title!")
    } else if (pct_in_title) {
      pct <- round((length(degs)/nrow(results))*100, 2)
      pct <- paste0(pct, " %")
      deg_list[[paste0(set_name, ": ", pct, " DE")]] <- degs
    } else if (n_in_title) {
      deg_list[[paste0(set_name, ": ", length(degs), " DEGs")]] <- degs
    } else {
      deg_list[[set_name]] <- degs
    }
  }
  
  return(deg_list)
}

# Generate set list for upset with basic set names (set pct_in_title or n_in_title to TRUE to include percent DE or number of DEGs in set names)
deg_list <- load_degs(result_dir)

# Generate initial plot to find total intersections
plt <- UpSetR::upset(UpSetR::fromList(deg_list), nsets = length(deg_list))
n_bars <- nrow(unique(plt$New_data))

# Generate colors for all intersections - adjust seed until colors look good (use any preferred method to generate bar colors)
set.seed(125)
bar_colors <- randomcoloR::distinctColorPalette(n_bars)

# Generate basic plot using UpSetR
plt <- UpSetR::upset(fromList(deg_list), order.by = "freq", nsets = length(deg_list), nintersects = 35,
                     sets.bar.color = "darkblue", main.bar.color = bar_colors)

# Optionally, generate plot using ComplexUpset (more customization options)

# Define function to get data for ComplexUpset (https://stackoverflow.com/questions/75990226/add-labels-to-upset-plot-so-the-values-of-intersection-would-be-visible-along-t)
from_list <- function(list_data) {
  members = unique(unlist(list_data))
  data.frame(
    lapply(list_data, function(set) members %in% set),
    row.names=members,
    check.names=FALSE
  )
}

# Define number of bars to show
n_bars <- 12

# Generate colors only for specified number of bars - adjust seed until colors look good (use any preferred method to generate bar colors)
set.seed(125)
bar_colors <- randomcoloR::distinctColorPalette(n_bars)

# Generate plot
plt <- ComplexUpset::upset(matrix,
                           intersect = colnames(matrix)[colnames(matrix) != "gene_name"], # Use all sets
                           base_annotations=list(
                             'Intersection size'= (
                               intersection_size(
                                 bar_number_threshold = 1,
                                 fill = bar_colors, # Colors will not show in exact order of bar_colors, but will be distinct
                                 width = 0.7,  # Adjust bar width
                                 text = list(size = 4.5) # Adjust size of labels above bars
                               ) 
                             )
                           ),
                           keep_empty_groups = TRUE, # Show all sets even if not included in any displayed intersections https://krassowski.github.io/complex-upset/reference/upset_data.html
                           matrix = intersection_matrix(geom = geom_point(size = 3.5)), # Adjust point size 
                           width_ratio = 0.3, # Increase to make set size bars wider relative to bar plot 
                           n_intersections = n_bars,
                           height_ratio = 0.5, # Increase to make point matrix taller relative to bar plot
                           set_sizes = (upset_set_size(filter_intersections = FALSE, geom = geom_bar(fill = "darkblue"))), # Color of set size bars
                           sort_sets = "ascending", # Sort sets in ascending order 
                           name = "Group", # General description of groups
                           themes = upset_default_themes(text = element_text(size = 12),# Adjust global text size (except for intersection size labels)
                                                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                         panel.background = element_blank())
)



