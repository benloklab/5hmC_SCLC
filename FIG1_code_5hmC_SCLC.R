##############################################################################
#
# Code found in this script was used to make figure 1 in the publication:
# "Exploring the Utility of Cell-Free DNA Hydroxymethylation Profiling in 
# Small-Cell Lung Cancer"
#
# This code can be used for creating boxplots, gene feature pie charts, 
# heatmaps, PCA, volcano plots, and dot plots to visualize GSEA data
# 
# Code was created by: Janice Li
# Last updated: May 14, 2026
#
##############################################################################


### Load libraries ---------------------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
library(ggpubr)

# Fig 1E
library(pheatmap)



### Load data --------------------------------------------------------------
# Data used in this manuscript can be found at: https://doi.org/10.5281/zenodo.19216411

# For counts matrices:
# * Rows = genomic regions/5hmC peaks
# * Columns = sample

# SCLC cfDNA 5hmC data and corresponding CDX gDNA 5hmC counts matrix
sclc_cdx_5hmC_cm <- readRDS("hmC.SCLC_CDX.count_mat.rds")

# SCLC and NCC cfDNA 5hmC counts matrix
sclc_ncc_5hmC_cm <- readRDS("hmC.SCLC_NCC.count_mat.rds")

# RNA-seq counts matrix (normalized to CPM)
rna_cdx <- readRDS("RNA.CDX.CPM_matrix.rds")

# NCC metadata
ncc_annot <- read.csv("NCC_clinical_annotations.csv")

# SCLC metadata
sclc_annot <- read.csv("SCLC_clinical_annotations.csv")


### Fig 1A, B, E ---------------------------------------------------------------

# Fig 1A was created in biorender.com
# Fig 1B was created in Microsoft Excel
# Fig 1E was created in dhmr.R code. Please go there for the full code for this 
  # figure. Note: Some input files may be very large and memory-intensive. We 
  # recommend running the analysis in a cloud or HPC environment rather than 
  # locally.


### Fig 1C - Global 5hmC boxplot --------------------------------------------

### [PREP] Variables for saving ###

# For saving RPKM medians
OUTPUT_PATH_med <- "<INSERT PATH/>"
FILE_NAME_med <- "<INSERT FILE_NAME>"

# For saving final boxplot
OUTPUT_PATH_p <- "<INSERT PATH/>"
FILE_NAME_p <- "<INSERT FILE_NAME>"

# Box plot colours
colours <- c("SCLC"="#DE5057", "NCC"="#4568BF") 


### [CODE] ###

### 1) Calculate median counts for each sample, across all 5hmC peak regions
rpkm_med <- apply(sclc_ncc_5hmC_cm, 2, median, na.rm = TRUE)

# [OPTIONAL] Save RPKM medians
OUTPUT_FILE_NAME_med <- paste0(OUTPUT_PATH_med, FILE_NAME_med)
saveRDS(rpkm_med, OUTPUT_FILE_NAME_med)


### 2) Add global median data to metadata

# Make new data frame with sample names
df <- data.frame("Samples" = names(rpkm_med), "rpkm_med" = rpkm_med)

# Get sample type
df$Type <- sub("\\..*", "", colnames(df))
df$Type <- factor(df$Type, levels = c("SCLC", "NCC"))

# Make new df with col medians, sample name, and sample type
new_samples_df <- samples_df[, c("sample_names", "Type", "med_rpkm")]


### 3) Plot Box plot
p <- ggplot(df, aes(x=Type, y=rpkm_med)) +
  geom_boxplot(aes(fill=Type), alpha=0.6, outliers=FALSE) + #geom_violin for violin
  
  # Points
  geom_point(aes(fill=Type), size=2, shape=21, 
             position=position_jitterdodge(0.2)) +
  
  # Modify titles
  labs(x = "Sample Type",
       y = "Median 5hmC RPKM") +
  
  # Change aesthetics (i.e. background theme, colour)
  theme(legend.position="none") +
  scale_fill_manual(values=colours, name="Sample Type") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        legend.position="none",
        text = element_text(size = 12),           # General text size
        axis.title = element_text(size = 14),     # Axis titles
        axis.text = element_text(size = 12)) +    # Axis tick labels
  
  # Calculate p-value (T-test for normally distributed data; 
  # wilcox test for non-normally distributed data)
  stat_compare_means(comparisons = list(c("SCLC", "NCC")),
                     method = "wilcox.test", label = "p.value", 
                     textsize = 8, vjust = 0.2) +
  
  scale_y_continuous(breaks = seq(0, max(df$med_rpkm, na.rm = TRUE), by = 0.5))

# Display plot
p

# Save plot
OUTPUT_FILE_NAME_p <- paste0(OUTPUT_PATH_p, FILE_NAME_p)
ggsave(OUTPUT_FILE_NAME_p, plot=p, width=3, height=4, dpi=300, units="in")



### [BONUS] Get median values
med_NCC <- median(df$rpkm_med[df$Type == "NCC"], na.rm = TRUE)
med_SCLC <- median(df$rpkm_med[df$Type == "SCLC"], na.rm = TRUE)

cat(paste0("The median of median 5hmC RPKMs for NCC is: ", round(med_NCC, 2), "\n"))
cat(paste0("The median of median 5hmC RPKMs for SCLC is: ", round(med_SCLC, 2), "\n"))

# cat("NCC\n")
# round(summary(df$med_rpkm[df$Type == "NCC"]), 2)
# 
# cat("SCLC\n")
# round(summary(df$med_rpkm[df$Type == "SCLC"]), 2)



### Fig 1D - Gene ft plots ----------------------------------------------------

### Please run annotation code "annot.R" before running this code! ###
# This will give you a file called "basic_genes_counts.csv"

### [PREP] ###

# Plot colours
colours <-  c( "#28405c", "#42899b","#94c4c1","#edada3", "#ce6a6c", "#cc222b")

# For saving final plot
OUTPUT_PATH_p <- "<INSERT PATH/>"
FILE_NAME_p <- "<INSERT FILE_NAME>"


### [CODE] ###

# 1) Import annotated basic gene feature annotations
gene_ft_counts <- read.csv("basic_genes_counts.csv", header=TRUE)

# Remove first column
gene_ft_counts <- gene_ft_counts[, -1]

# Re-label
gene_ft_counts$type <- sapply(strsplit(gene_ft_counts$sample, "\\."), function(x) x[2])
gene_ft_counts$feature <- factor(gene_ft_counts$feature, 
                                 levels = c("hg19_genes_introns",
                                            "hg19_genes_exons",
                                            "hg19_genes_1to5kb", 
                                            "hg19_genes_promoters",
                                            "hg19_genes_3UTRs",
                                            "hg19_genes_5UTRs"))
gene_ft_counts$type <- factor(gene_ft_counts$type, 
                              levels = c("SCLC", "NCC"))

gene_ft <- c("Introns","Exons","1 to 5 kb", "Promoters", "3\' UTR", "5\' UTR")


# 2) Calculate average counts for each gene feature
med_counts <- gene_ft_counts %>%
  group_by(feature, type) %>%
  summarize(med_counts = median(count)) %>%
  ungroup()


# Calculate total counts by Type
total_by_type <- tapply(med_counts$med_counts, med_counts$type, sum)


# Add a new column with percentage by Type
med_counts$percent <- round(med_counts$med_counts / total_by_type[med_counts$type] * 100, 0)


### [PLOT] 3) Plot stacked bar plot ###
p <- ggplot(med_counts, aes(x=type, y=med_counts, fill=feature, group=type)) +  
  
  geom_bar(position="fill", 
           stat="identity", 
           color="black",
           linewidth=0.1, 
           alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_manual(values=colours, labels = gene_ft) +
  ggtitle("Gene feature plot") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_blank()) +
  guides(fill=guide_legend(reverse=TRUE))+
  xlab("Sample Type") +
  ylab("Percent") +
  geom_text(aes(label = paste0(percent, " %")),
            position = position_fill(vjust = 0.5),
            size = 5,
            color = "black")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))  

# Visualize plot
p

# Save plot
FILE_PATH <- paste0(OUTPUT_PATH_p, FILE_NAME_p)
ggsave(FILE_PATH, plot=p, width=3.5, height=5, dpi=300, units="in")


### 3) [PLOT] Plot pie chart ###
med_counts <- med_counts[med_counts$type == "SCLC",] # Run this for SCLC
# med_counts <- med_counts[med_counts$type == "NCC",] # Run this for NCC


# Plotting gene feature plot 
p <- ggplot(med_counts_sclc, aes(x=type, y=med_counts, fill=feature, group=type)) +  
  geom_bar(position="fill", 
           stat="identity", 
           color="black",
           linewidth=0.1, 
           alpha=0.7) +
  coord_polar("y", start = 0) +
  theme(legend.position="none") +
  scale_fill_manual(values=colours, labels = gene_ft) +
  ggtitle("Gene feature plot") +
  theme_void() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 10, l = 5, unit = "mm"),  
        legend.title = element_blank(),
        legend.position = "bottom") +
  guides(fill=guide_legend(reverse=TRUE)) +
  geom_text(aes(label = paste0(percent, " %")),
            position = position_fill(vjust = 0.5),
            size = 5,
            color = "black")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))  

# Visualize plot
p

# Save plot
FILE_PATH <- paste0(OUTPUT_PATH_p, FILE_NAME_p)
ggsave(FILE_PATH, plot=p, width=3.5, height=5, dpi=300, units="in")



### Fig 1F - PCA --------------------------------------------

### Please run annotation code "dhmr.R" before running this code! ###
# This will give you a file called "pca_data.rds"

### [PREP] ###

# Annotations
sclc_annot <- sclc_annot[, 1:4] # Subset to variables of interest
annot <- rbind(sclc_annot, ncc_annot)

# Load PCA data
pca_data <- readRDS("pca_data.rds")

# Setting colours
namedColor <- c("SCLC"="#E8858A", "NCC"="#7C96D2") 

# Variance percent of each principal component
percentVar <- round(100*attr(pca_data, "percentVar"))

# For saving final plot
OUTPUT_PATH_p <- "<INSERT PATH/>"
FILE_NAME_p <- "<INSERT FILE_NAME>"


### [CODE] ###

### Plot PCA
p <- ggplot(pca_data, aes(x=PC1, y=PC2)) +
  geom_point(size=4, shape=21, color="black", aes(fill=Type)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  scale_fill_manual(values=namedColor) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(fill="Sample Type") 

# Visualize PCA
p

# Save
FILE_PATH <- paste0(OUTPUT_PATH_p, FILE_NAME_p)
ggsave(FILE_PATH, plot=p, width=6, height=3, dpi=300, units="in")


### Fig 1G - Volcano ----------------------------------------
### Please run annotation code "dhmr.R" before running this code! ###
# This will give you a file called "res.rds", which are results after 
  # differential analysis

### [PREP] ###

# Import res and convert to a data frame
res <- readRDS("res.rds")
res_df <- data.frame(res)

## Visualize res_df
# dim(res_df) 
# head(res_df)


# Setting cut-offs
res_df$padj <- -log10(res_df$padj) # Transforms padj
pval.cutoff <- 0.05
fc.cutoff <- 1
res_df$colors <- "Nonsignificant"


# Setting volcano plot colours
namedColor <- c("NCC"="#7C96D2", "SCLC"="#E8858A", "Nonsignificant"="#D9D9D9")

# RED
res_df$colors[which(res_df$padj > -log10(pval.cutoff) & res_df$log2FoldChange > fc.cutoff)] <- "SCLC"
# length(which(res_df$colors == "SCLC"))

# BLUE
res_df$colors[which(res_df$padj > -log10(pval.cutoff) & res_df$log2FoldChange < -fc.cutoff)] <- "NCC"
# length(which(res_df$colors == "NCC"))


# For saving final plot
OUTPUT_PATH_p <- "<INSERT PATH/>"
FILE_NAME_p <- "<INSERT FILE_NAME>" # RECOMMENDED to save as a jpg and not a vector-based file!


### [CODE] ###

### 1) Clean res_df by taking out NA/inf log2FC
plot_df <- res_df %>%
  filter(!is.na(-log10(res_df$padj)), is.finite(log2FoldChange))

### 2) Plot
p <- ggplot(plot_df, aes(x=log2FoldChange, y=padj, color=colors)) +
  geom_point(shape=19) +
  geom_hline(yintercept = -log10(pval.cutoff), color="black", lty=2, lwd=0.5) +
  geom_vline(xintercept = c(-fc.cutoff, fc.cutoff), color="black", lty=2, lwd=0.5) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  scale_color_manual(values=namedColor) +
  labs(x="log2(fold change)",
       y="-log10(adjusted p-value)")

# Visualize plot
p

# Save plot
FILE_PATH <- paste0(OUTPUT_PATH_p, FILE_NAME_p)
ggsave(FILE_PATH, plot=p, width=6, height=3, dpi=300, units="in")



### Fig 1H - GSEA -------------------------------------------
### Please run annotation code "gsea.R" before running this code! ###
# This will give you a file called "fgsea_res.csv", which are results after 
# differential analysis

### [PREP] ###
# Import fgsea data
fgsea_go_bp <- read.csv("fgsea_res.csv")

# Define colors
nes_colors <- c(
  "FALSE" = "#AED4DD",   # negative NES
  "TRUE"  = "#E9AFC6"    # positive NES
)

# For saving final plot
OUTPUT_PATH_p <- "<INSERT PATH/>"
FILE_NAME_p <- "<INSERT FILE_NAME>"


### [CODE] ###

### 1) Reorganize fgsea resuts and select top 10 pathways

# Order fgsea results
fgsea_go_bp <- fgsea_go_bp[order(fgsea_go_bp$padj, -abs(fgsea_go_bp$NES)), ]

# Select top 10 pathways
top_pathways <- fgsea_go_bp[1:10, ]

# Clean pathway names
clean_terms <- gsub("_", " ", gsub("^GOBP_", "", top_pathways$pathway))
top_pathways$short <- clean_terms
top_pathways$short <- stringr::str_wrap(top_pathways$short, width = 30)

# Convert to factor so order is preserved
top_pathways$short <- factor(
  top_pathways$short,
  levels = top_pathways$short[order(top_pathways$NES)]
)


### 2) Plot dot plot
p <- ggplot(
  top_pathways,
  aes(
    x = short,
    y = NES,
    fill = NES > 0
  )
) +
  
  geom_point(
    shape = 21,           
    size = 5,
    color = "black",     
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
    axis.title.y = element_text(margin = margin(r = 10))
  ) +
  
  guides(fill = "none")

# Visualize
p

# Save
FILE_PATH <- paste0(OUTPUT_PATH_p, FILE_NAME_p)
ggsave(FILE_PATH, plot=p, width=6, height=5, dpi=300, units="in")
