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
# Summary: Compare average QC metrics per donor per brain region for lecanemab vs. nAD using Wilcoxon test
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
  library("ggplot2")
  library("randomcoloR")
})

# Define output folder
output_dir <- "/path/to/qc/barplots/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load integrated cohort 5/7/8 Seurat object
s <- readRDS("/path/to/integrated/cohort578/object.rds")

# Define condition groups
s$group <- NA
s$group[grep("\\.B", s$sample_id)] <- "LCMB"
s$group[grep("\\.A|^A", s$sample_id)] <- "CAA"

# Define filter operator 
`%notin%` <- Negate(`%in%`)

# Define variable for brain region
s@meta.data <- s@meta.data %>% mutate(region = case_when(substr(s$sample_id, nchar(s$sample_id), nchar(s$sample_id)) == 1 ~ "FCX",
                                                         substr(s$sample_id, nchar(s$sample_id), nchar(s$sample_id)) == 3 ~ "TCX",
                                                         substr(s$sample_id, nchar(s$sample_id), nchar(s$sample_id)) == 4 ~ "PCX",
                                                         substr(s$sample_id, nchar(s$sample_id), nchar(s$sample_id)) == 9 ~ "HIPP"))

# Define broad donor variable
s$donor <- NA
s$donor[grep("^A14.193", s$sample_id)] <- "A14.193"
s$donor[grep("^A11.170", s$sample_id)] <- "A11.170"
s$donor[grep("^NMA22.A", s$sample_id)] <- "NMA22.A"
s$donor[grep("^NMA22.B", s$sample_id)] <- "NMA22.B"
print(unique(s@meta.data[,c("group", "region", "sample_id", "donor")]))

# Extract meta data
data <- s@meta.data

# Calculate average MT% and nFeatures per donor per region 
data <- data %>% dplyr::group_by(group, donor, region) %>% dplyr::summarize(mt = mean(percent.mt), nfeat = mean(nFeature_Spatial))

# Set factor levels for region and donor
data$region <- factor(data$region, levels = c("FCX", "TCX", "PCX", "HIPP"))
data$donor <- factor(data$donor, levels = c("A14.193", "A11.170", "NMA22.A", "NMA22.B"))

# Define color for region
colors <- c("#8DD3C7", "#FDB462", "#3399CC", "#E78AC3")
names(colors) <- c("FCX", "TCX", "PCX", "HIPP")

# MT %

# Calculate mean and SEM per group 
bar_data <- data %>% dplyr::group_by(donor) %>% dplyr::summarize(mean = mean(mt), sem = plotrix::std.error(mt))
bar_data$donor <- factor(bar_data$donor, levels = c("A14.193", "A11.170", "NMA22.A", "NMA22.B"))
bar_data$lower <- bar_data$mean - bar_data$sem
bar_data$upper <- bar_data$mean + bar_data$sem

plt <- ggplot(data, aes(x = donor, y = mt, color = region)) + 
  geom_col(data = bar_data, mapping = aes(x = donor, y = mean), inherit.aes = FALSE, 
           alpha = 0.5, width = 0.5, size = 1.5, color = "#505050", fill = "lightgray") +
  geom_errorbar(data = bar_data, mapping = aes(x = donor, ymin = lower, ymax = upper), inherit.aes = FALSE, 
                width = 0.1, size = 1.5, color = "#505050") +
  geom_jitter(stat="identity", height = 0, width = 0.1, size = 3) +
  labs(color = NULL, y = "Percent MT Expression", x = NULL) + theme_classic() + 
  theme(plot.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, max(data$mt)+0.5)) +
  scale_color_manual(values = colors) +
  theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pdf(paste0(output_dir, "mt_barplot.pdf"), height = 8, width = 8)
print(plt)
dev.off()

# Wilcoxon test per region (exact p values used to define significance)
stats <- data.frame() 
for (region in levels(data$region)) {
  cur_data <- data[data$region == region,]
  exact <- stats::wilcox.test(x = cur_data$mt[cur_data$group == "LCMB"], y = cur_data$mt[cur_data$group == "CAA"],
                              alternative = "two.sided", paired = FALSE, exact = TRUE, correct = FALSE)
  adj <- stats::wilcox.test(x = cur_data$mt[cur_data$group == "LCMB"], y = cur_data$mt[cur_data$group == "CAA"],
                            alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE)
  print(exact$statistic == adj$statistic)
  row <- data.frame(region = region, w_stat = exact$statistic, p_exact = exact$p.value, p_corrected = adj$p.value)
  stats <- rbind(stats, row)
}
write.csv(stats, paste0(output_dir, "mt_wilcoxon.csv"), row.names = FALSE)

# nFeatures

# Calculate mean and SEM per group 
bar_data <- data %>% dplyr::group_by(donor) %>% dplyr::summarize(mean = mean(nfeat), sem = plotrix::std.error(nfeat))
bar_data$donor <- factor(bar_data$donor, levels = c("A14.193", "A11.170", "NMA22.A", "NMA22.B"))
bar_data$lower <- bar_data$mean - bar_data$sem
bar_data$upper <- bar_data$mean + bar_data$sem

plt <- ggplot(data, aes(x = donor, y = nfeat, color = region)) + 
  geom_col(data = bar_data, mapping = aes(x = donor, y = mean), inherit.aes = FALSE, 
           alpha = 0.5, width = 0.5, size = 1.5, color = "#505050", fill = "lightgray") +
  geom_errorbar(data = bar_data, mapping = aes(x = donor, ymin = lower, ymax = upper), inherit.aes = FALSE, 
                width = 0.1, size = 1.5, color = "#505050") +
  geom_jitter(stat="identity", height = 0, width = 0.1, size = 3) +
  labs(color = NULL, y = "Genes Detected", x = NULL) + theme_classic() + 
  theme(plot.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, max(data$nfeat)+30)) +
  scale_color_manual(values = colors) +
  theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pdf(paste0(output_dir, "nfeat_barplot.pdf"), height = 8, width = 8)
print(plt)
dev.off()

# Wilcoxon test per region (exact p values used to define significance)
stats <- data.frame() 
for (region in levels(data$region)) {
  cur_data <- data[data$region == region,]
  exact <- stats::wilcox.test(x = cur_data$nfeat[cur_data$group == "LCMB"], y = cur_data$nfeat[cur_data$group == "CAA"],
                              alternative = "two.sided", paired = FALSE, exact = TRUE, correct = FALSE)
  adj <- stats::wilcox.test(x = cur_data$nfeat[cur_data$group == "LCMB"], y = cur_data$nfeat[cur_data$group == "CAA"],
                            alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE)
  print(exact$statistic == adj$statistic)
  row <- data.frame(region = region, w_stat = exact$statistic, p_exact = exact$p.value, p_corrected = adj$p.value)
  stats <- rbind(stats, row)
}
write.csv(stats, paste0(output_dir, "nfeat_wilcoxon.csv"), row.names = FALSE)
