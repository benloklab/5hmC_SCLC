##############################################################################
# Figure 1 Visualization Code
#
# Publication:
#   Exploring the Utility of Cell-Free DNA Hydroxymethylation Profiling in
#   Small-Cell Lung Cancer
#
# Description:
#   Code for generating Figure 1 visualizations, including global 5hmC boxplots,
#   gene feature plots, PCA, volcano plots, and GSEA dot plots.
#
# Original author: Janice Li
# Last updated: May 14, 2026
##############################################################################

### Libraries --------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(stringr)

### Input Data -------------------------------------------------------------------
# Data used in this manuscript can be found at:
# https://doi.org/10.5281/zenodo.19216411
#
# Counts matrix format:
#   Rows    = genomic regions / 5hmC peaks
#   Columns = samples

# Counts matrices
sclc_cdx_5hmC_cm <- readRDS("hmC.SCLC_CDX.count_mat.rds")
sclc_ncc_5hmC_cm <- readRDS("hmC.SCLC_NCC.count_mat.rds")
rna_cdx          <- readRDS("RNA.CDX.CPM_matrix.rds")

# Annotation data
ncc_annot  <- read.csv("NCC_clinical_annotations.csv")
sclc_annot <- read.csv("SCLC_clinical_annotations.csv")

# Figure Notes -----------------------------------------------------------------
# Fig. 1A was created in BioRender.
# Fig. 1B was created in Microsoft Excel.
# Fig. 1E was created in dhmr.R. See that script for the full code.
#
# Note:
#   Some input files may be very large and memory-intensive. Running the analysis
#   in a cloud or HPC environment is recommended.

### Figure 1C - Global 5hmC Boxplot ----------------------------------------------

# Output variables
output_path_medians <- "<INSERT PATH/>"
output_file_medians <- "<INSERT FILE_NAME>"

output_path_plot <- "<INSERT PATH/>"
output_file_plot <- "<INSERT FILE_NAME>"

# Plot colors
boxplot_colors <- c(
  "SCLC" = "#DE5057",
  "NCC"  = "#4568BF"
)

### 1) Calculate median counts for each sample, across all 5hmC peak regions
rpkm_medians <- apply(sclc_ncc_5hmC_cm, 2, median, na.rm = TRUE)

# [OPTIONAL] Save RPKM medians
saveRDS(
  rpkm_medians,
  file = file.path(output_path_medians, output_file_medians)
)

### 2) Add global median data to metadata
# Make new data frame with sample names
global_5hmc_df <- data.frame(
  sample   = names(rpkm_medians),
  rpkm_med = rpkm_medians
)

# Get sample type
global_5hmc_df$Type <- sub("\\..*", "", global_5hmc_df$sample)
global_5hmc_df$Type <- factor(global_5hmc_df$Type, levels = c("SCLC", "NCC"))

### 3) Plot Box plot
fig_1c <- ggplot(global_5hmc_df, aes(x = Type, y = rpkm_med)) +
  geom_boxplot(
    aes(fill = Type),
    alpha = 0.6,
    outliers = FALSE
  ) +
  
  # Points
  geom_point(
    aes(fill = Type),
    size = 2,
    shape = 21,
    position = position_jitterdodge(jitter.width = 0.2)
  ) +
  
  # Calculate p-value
  stat_compare_means(
    comparisons = list(c("SCLC", "NCC")),
    method = "wilcox.test",
    label = "p.value",
    textsize = 8,
    vjust = 0.2
  ) +
  
  # Change aesthetics
  scale_fill_manual(values = boxplot_colors, name = "Sample Type") +
  scale_y_continuous(
    breaks = seq(
      0,
      max(global_5hmc_df$rpkm_med, na.rm = TRUE),
      by = 0.5
    )
  ) +
  labs(
    x = "Sample Type",
    y = "Median 5hmC RPKM"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    text             = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12)
  )

# Display plot
fig_1c

# Save plot
ggsave(
  filename = file.path(output_path_plot, output_file_plot),
  plot     = fig_1c,
  width    = 3,
  height   = 4,
  dpi      = 300,
  units    = "in"
)

# Summary of medians
median_ncc  <- median(global_5hmc_df$rpkm_med[global_5hmc_df$Type == "NCC"], na.rm = TRUE)
median_sclc <- median(global_5hmc_df$rpkm_med[global_5hmc_df$Type == "SCLC"], na.rm = TRUE)

cat("Median of median 5hmC RPKMs for NCC: ", round(median_ncc, 2), "\n")
cat("Median of median 5hmC RPKMs for SCLC:", round(median_sclc, 2), "\n")

### Fig 1D - Gene feature Plots ----------------------------------------------------
# Before running:
#   Run annot.R first. This should create "basic_genes_counts.csv".

# Plot colours
gene_feature_colors <- c(
  "#28405c",
  "#42899b",
  "#94c4c1",
  "#edada3",
  "#ce6a6c",
  "#cc222b"
)

# Gene feature labels
gene_feature_labels <- c(
  "Introns",
  "Exons",
  "1 to 5 kb",
  "Promoters",
  "3' UTR",
  "5' UTR"
)

# Output variables
output_path_plot <- "<INSERT PATH/>"
output_file_plot <- "<INSERT FILE_NAME>"

### 1) Import annotated basic gene feature annotations
gene_feature_counts <- read.csv("basic_genes_counts.csv", header = TRUE)

# Remove first column if it is an index column
gene_feature_counts <- gene_feature_counts[, -1]

### 2) Clean labels
gene_feature_counts$type <- sapply(
  strsplit(gene_feature_counts$sample, "\\."),
  function(x) x[2]
)

gene_feature_counts$feature <- factor(
  gene_feature_counts$feature,
  levels = c(
    "hg19_genes_introns",
    "hg19_genes_exons",
    "hg19_genes_1to5kb",
    "hg19_genes_promoters",
    "hg19_genes_3UTRs",
    "hg19_genes_5UTRs"
  )
)

gene_feature_counts$type <- factor(
  gene_feature_counts$type,
  levels = c("SCLC", "NCC")
)

### 3) Calculate median counts and percentages
median_gene_feature_counts <- gene_feature_counts %>%
  group_by(feature, type) %>%
  summarize(med_counts = median(count, na.rm = TRUE), .groups = "drop")

total_by_type <- tapply(
  median_gene_feature_counts$med_counts,
  median_gene_feature_counts$type,
  sum
)

median_gene_feature_counts$percent <- round(
  median_gene_feature_counts$med_counts /
    total_by_type[median_gene_feature_counts$type] * 100,
  0
)

### 4) Pie chart
# Change selected_type to "NCC" for the NCC pie chart.
selected_type <- "SCLC"

median_gene_feature_counts_pie <- median_gene_feature_counts %>%
  filter(type == selected_type)

fig_1d_pie <- ggplot(
  median_gene_feature_counts_pie,
  aes(x = type, y = med_counts, fill = feature, group = type)
) +
  geom_bar(
    position = "fill",
    stat     = "identity",
    color    = "black",
    linewidth = 0.1,
    alpha    = 0.7
  ) +
  coord_polar("y", start = 0) +
  geom_text(
    aes(label = paste0(percent, " %")),
    position = position_fill(vjust = 0.5),
    size     = 5,
    color    = "black"
  ) +
  scale_fill_manual(values = gene_feature_colors, labels = gene_feature_labels) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(title = paste("Gene Feature Plot:", selected_type)) +
  theme_void() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(t = 5, r = 5, b = 10, l = 5, unit = "mm"),
    legend.title     = element_blank(),
    legend.position  = "bottom"
  )

# Visualize plot
fig_1d_pie

# Save plot
ggsave(
  filename = file.path(output_path_plot, output_file_plot),
  plot     = fig_1d_pie,
  width    = 3.5,
  height   = 5,
  dpi      = 300,
  units    = "in"
)

### [OPTIONAL: Stacked bar plot]
fig_1d_bar <- ggplot(
  median_gene_feature_counts,
  aes(x = type, y = med_counts, fill = feature, group = type)
) +
  geom_bar(
    position = "fill",
    stat     = "identity",
    color    = "black",
    linewidth = 0.1,
    alpha    = 0.7
  ) +
  geom_text(
    aes(label = paste0(percent, " %")),
    position = position_fill(vjust = 0.5),
    size     = 5,
    color    = "black"
  ) +
  scale_fill_manual(values = gene_feature_colors, labels = gene_feature_labels) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(
    title = "Gene Feature Plot",
    x     = "Sample Type",
    y     = "Percent"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title     = element_blank(),
    legend.position  = "none"
  )

# Visualize plot
fig_1d_bar

# Save plot
ggsave(
  filename = file.path(output_path_plot, output_file_plot),
  plot     = fig_1d_bar,
  width    = 3.5,
  height   = 5,
  dpi      = 300,
  units    = "in"
)

### Fig 1F - PCA --------------------------------------------
# Before running:
#   Run dhmr.R first. This should create "pca_data.rds".

# Annotations
sclc_annot_subset <- sclc_annot[, 1:4]
sample_annotations <- rbind(sclc_annot_subset, ncc_annot)

# Load PCA data
pca_data <- readRDS("pca_data.rds")

# Setting colours
pca_colors <- c(
  "SCLC" = "#E8858A",
  "NCC"  = "#7C96D2"
)

# Percent variance for each principal component
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Output variables
output_path_plot <- "<INSERT PATH/>"
output_file_plot <- "<INSERT FILE_NAME>"

### 1) Plot PCA
fig_1f <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 4, shape = 21, color = "black", aes(fill = Type)) +
  scale_fill_manual(values = pca_colors) +
  labs(
    x    = paste0("PC1: ", percent_var[1], "% variance"),
    y    = paste0("PC2: ", percent_var[2], "% variance"),
    fill = "Sample Type"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Visualize PCA
fig_1f

# Save plot
ggsave(
  filename = file.path(output_path_plot, output_file_plot),
  plot     = fig_1f,
  width    = 6,
  height   = 3,
  dpi      = 300,
  units    = "in"
)

### Fig 1G - Volcano ----------------------------------------
# Before running:
#   Run dhmr.R first. This should create "res.rds".

res <- readRDS("res.rds")
res_df <- data.frame(res)

# Set cut-offs
pval_cutoff <- 0.05
fc_cutoff   <- 1

# Setting volcano plot colours
volcano_colors <- c(
  "NCC"            = "#7C96D2",
  "SCLC"           = "#E8858A",
  "Nonsignificant" = "#D9D9D9"
)

res_df <- res_df %>%
  mutate(
    neg_log10_padj = -log10(padj),
    colors = case_when(
      neg_log10_padj > -log10(pval_cutoff) & log2FoldChange >  fc_cutoff ~ "SCLC",
      neg_log10_padj > -log10(pval_cutoff) & log2FoldChange < -fc_cutoff ~ "NCC",
      TRUE ~ "Nonsignificant"
    )
  ) %>%
  filter(!is.na(neg_log10_padj), is.finite(log2FoldChange))

# Output variables
output_path_plot <- "<INSERT PATH/>"
output_file_plot <- "<INSERT FILE_NAME>" # Recommended: save as JPG, not vector.

### 1) Plot volcano
fig_1g <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = colors)) +
  geom_point(shape = 19) +
  geom_hline(
    yintercept = -log10(pval_cutoff),
    color      = "black",
    lty        = 2,
    lwd        = 0.5
  ) +
  geom_vline(
    xintercept = c(-fc_cutoff, fc_cutoff),
    color      = "black",
    lty        = 2,
    lwd        = 0.5
  ) +
  scale_color_manual(values = volcano_colors) +
  labs(
    x = "log2(fold change)",
    y = "-log10(adjusted p-value)"
  ) +
  theme_bw() +
  theme(
    legend.title     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Visualize plot
fig_1g

# Save plot
ggsave(
  filename = file.path(output_path_plot, output_file_plot),
  plot     = fig_1g,
  width    = 6,
  height   = 3,
  dpi      = 300,
  units    = "in"
)

### Fig 1H - GSEA -------------------------------------------
# Before running:
#   Run gsea.R first. This should create "fgsea_res.csv".

# Import fgsea data
fgsea_go_bp <- read.csv("fgsea_res.csv")

# Define colors
nes_colors <- c(
  "FALSE" = "#AED4DD",
  "TRUE"  = "#E9AFC6"
)

# Output variables
output_path_plot <- "<INSERT PATH/>"
output_file_plot <- "<INSERT FILE_NAME>"

### 1) Reorganize fgsea resuts and select top 10 pathways
top_pathways <- fgsea_go_bp %>%
  arrange(padj, desc(abs(NES))) %>%
  slice_head(n = 10) %>%
  mutate(
    short = pathway %>%
      gsub(pattern = "^GOBP_", replacement = "") %>%
      gsub(pattern = "_", replacement = " ") %>%
      str_wrap(width = 30),
    short = factor(short, levels = short[order(NES)])
  )

### 2) Plot dot plot
fig_1h <- ggplot(top_pathways, aes(x = short, y = NES, fill = NES > 0)) +
  geom_point(
    shape  = 21,
    size   = 5,
    color  = "black",
    stroke = 0.6
  ) +
  scale_fill_manual(values = nes_colors) +
  coord_flip() +
  labs(
    x = "GO Biological Process",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y     = element_text(margin = margin(r = 10))
  ) +
  guides(fill = "none")

# Visualize plot
fig_1h

# Save plot
ggsave(
  filename = file.path(output_path_plot, output_file_plot),
  plot     = fig_1h,
  width    = 6,
  height   = 5,
  dpi      = 300,
  units    = "in"
)
