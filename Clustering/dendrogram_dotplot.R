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
# Summary: Add dendrogram to Seurat dot plot of marker gene expression
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
  library("cowplot")
})

# Load Seurat object 
s <- readRDS("/path/to/Seurat/object.rds")

# Recorrect SCT data
if ("Spatial" %in% names(s@assays)) {
  for (name in names(s@assays$SCT@SCTModel.list)) {
    s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
  }
}
DefaultAssay(s) <- "SCT"
s <- PrepSCTFindMarkers(s)

# Define cluster variable and set idents
s$cluster <- s@meta.data[["<name of cluster column>"]]
Idents(s) <- "cluster"

# Load FindMarkers output 
markers <- read.csv("/path/to/markers.csv")

# Define marker genes to plot
marker_genes_to_plot <- c("<list of marker genes>")

# Optionally, select top n markers per cluster based on probability fold change 
n <- 2
marker_genes_to_plot <- c()
for(clust in all_markers$cluster|>unique()){
  data <- dplyr::filter(all_markers, cluster==clust)
  print(sum(data$PFC == Inf))
  data <- data %>% dplyr::arrange(desc(PFC))
  marker_genes_to_plot <- c(marker_genes_to_plot, data$gene|>head(n))
}
marker_genes_to_plot <- marker_genes_to_plot %>% unique()

# Modified original dot plot function return value to include more elements https://github.com/satijalab/seurat/blob/HEAD/R/visualization.R
dot_plot <- function(
    object,
    features,
    assay = NULL,
    cols = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = 'radius',
    scale.min = NA,
    scale.max = NA
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, cells = colnames(object[[assay]]), idents = idents))
  data.features <- FetchData(object = object, vars = features, cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- FetchData(object = object, vars = split.by)[cells, split.by]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop(paste0("Need to specify at least ", length(x = unique(x = splits)), " colors using the cols parameter"))
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning(
      "Scaling data with a low number of groups may produce misleading results",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = log1p(data.use))
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- unlist(x = lapply(
      X = data.plot$id,
      FUN = function(x)
        sub(
          paste0(".*_(",
                 paste(sort(unique(x = splits), decreasing = TRUE),
                       collapse = '|'
                 ),")$"),
          "\\1",
          x
        )
    )
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  }
  color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = 'Features',
      y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
    ) +
    theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(
      facets = ~feature.groups,
      scales = "free_x",
      space = "free_x",
      switch = "y"
    ) + theme(
      panel.spacing = unit(x = 1, units = "lines"),
      strip.background = element_blank()
    )
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
  }
  return(list(plot = plot, data.plot = data.plot, color.by = color.by, scale.min = scale.min, scale.max = scale.max,
              dot.scale = dot.scale, cols = cols, split.colors = split.colors, feature.groups = feature.groups,
              split.by = split.by, scale.func = scale.func))
}

# Extract data used for dot plot
full_output <- dot_plot(object = s, features = marker_genes_to_plot, assay = "SCT")
data.plot <- full_output[['data.plot']]
scaled <- data.plot[,c("features.plot", "id", "avg.exp.scaled")] %>% pivot_wider(names_from = "id", values_from = "avg.exp.scaled")
orig_cols <- colnames(scaled)
scaled <- data.frame(scaled)
colnames(scaled) <- orig_cols
rownames(scaled) <- scaled$features.plot
scaled$features.plot <- NULL

# Cluster columns (genes)
col_clust <- hclust(dist(as.matrix(scaled)))
colorder <- col_clust$labels[col_clust$order]
data.plot$features.plot <- as.character(data.plot$features.plot)
data.plot$features.plot <- factor(data.plot$features.plot, levels = colorder)

# Cluster rows (clusters)
row_clust <- hclust(dist(as.matrix(t(scaled))))
roworder <- row_clust$labels[row_clust$order]
data.plot$id <- as.character(data.plot$id)
data.plot$id <- factor(data.plot$id, levels = roworder)

# Generate continuous x and y variables
data.plot$y <- NA
i <- 1
for (cluster in levels(data.plot$id)) {
  data.plot$y[data.plot$id == cluster] <- i
  i <- i + 1
}

data.plot$x <- NA
i <- 1
for (gene in levels(data.plot$features.plot)) {
  data.plot$x[data.plot$features.plot == gene] <- i
  i <- i + 1
}


# Remake plot using code from function
scale.func <- full_output[['scale.func']]
plot <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) + # Updated x and y
  geom_point(mapping = aes_string(size = 'pct.exp', color = full_output[["color.by"]])) +
  scale.func(range = c(0, full_output[["dot.scale"]]), limits = c(full_output[["scale.min"]], full_output[["scale.max"]])) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(size = guide_legend(title = 'Percent Expressed')) +
  labs(
    x = 'Features',
    y = ifelse(test = is.null(x = full_output[["split.by"]]), yes = 'Identity', no = 'Split Identity')
  ) +
  theme_cowplot()
if (!is.null(x = full_output[["feature.groups"]])) {
  plot <- plot + facet_grid(
    facets = ~full_output[["feature.groups"]],
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) + theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  )
}
if (full_output[["split.colors"]]) {
  plot <- plot + scale_color_identity()
} else if (length(x = full_output[["cols"]]) == 1) {
  plot <- plot + scale_color_distiller(palette = full_output[["cols"]])
} else {
  plot <- plot + scale_color_gradient(low = full_output[["cols"]][1], high = full_output[["cols"]][2])
}
if (!full_output[["split.colors"]]) {
  plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
}

# Add dendrograms
plot <- plot + ggdendroplot::geom_dendro(row_clust, pointing = "side", xlim=c(31.5, 32.5))
plot <- plot + ggdendroplot::geom_dendro(col_clust, pointing = "updown", ylim=c(16.5, 17.5))

# Adjust x and y limits so that dots fit 
plot <- plot + scale_x_continuous(expand = c(0.02,0), breaks = 1:length(colorder), labels = colorder) + 
  scale_y_continuous(expand = c(0.02,0), breaks = 1:length(roworder), labels = roworder)
plot <- plot + theme(axis.text.x = element_text(angle = 45, size = 8, vjust = 1, hjust=1), axis.text.y = element_text(size = 8))





