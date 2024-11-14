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
# Summary: Generate density UMAPs for lecanemab and AN1792 amyloid-rich spots in gray matter
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
  library("randomcoloR")
  library("ggpubr")
  library("gridExtra")
  library("ggplot2")
  library("ggnewscale")
  library("ggheatmap")
  library("radiant.data")
  library("ggtree")
  library("aplot")
  library("ggrepel")
  library("ggplotify")
  library("grid")
  library("ggdendroplot")
})

# Define output folder 
output_dir <- "/path/to/plot/folder/"

# Load all-cohorts gray matter amyloid-rich integrated Seurat object 
s <- readRDS("/path/to/integrated/all/cohorts/amyloid/object.rds")

# Define function to calculate density
# https://slowkow.com/notes/ggplot2-color-by-density/
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Generate separate plots for iAD LCMB, nAD LCMB, iAD AN1792, nAD AN1792 
for (group in unique(s$condition)) {
  temp_s <- subset(s, condition == group)
  data <- temp_s@reductions$harmonyumap@cell.embeddings %>% as.data.frame
  data$density <- get_density(data$harmonyumap_1, data$harmonyumap_2, n = 100)
  
  plt <- ggplot(data, aes(x = harmonyumap_1, y = harmonyumap_2)) + 
    geom_point(aes(fill = density), alpha = 0.5, shape = 21, stroke = NA, size = 4) +
    stat_density_2d(geom = "density_2d", linewidth = 0.5, color = "gray30") +
    theme_classic() + viridis::scale_fill_viridis(option = "magma") +
    theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(aspect.ratio = 1)
  
  pdf(paste0(output_dir, group, ".pdf"), height = 10, width = 10)
  print(plt)
  dev.off()
}

