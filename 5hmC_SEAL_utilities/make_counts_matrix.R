##############################################################################
#
# Create a merged consensus peak BED file using DiffBind
#
# # This script identifies hydroxymethylation peak regions from MACS2 peak files
# and creates a merged consensus peak set for downstream featureCounts analysis.
#
# INPUT: csv file with column headers: Sample_ID, Condition, Replicate, bamReads, 
#   Peak, and Peakcaller
#
# OUTPUT: BED file containing merged consensus peaks
# 
# Reference:
#   https://hbctraining.github.io/Intro-to-ChIPseq/lessons/08_diffbind_differential_peaks.html
# 

# Author: Janice Li
# Last updated: May 21, 2026
#
##############################################################################

### Load libraries -------------------------------------------------------------
library(DiffBind)


### Load data ------------------------------------------------------------------


### User-defined paths ---------------------------------------------------------
sample_sheet <- "<PATH_TO_SAMPLE_SHEET.csv>"
  # Note that the sample sheet has 6 columns with the headers: 
    # Sample_ID: sample name
    # Condition: condition (i.e. cancer vs non-cancer)
    # Replicate: replicate number
    # bamReads: path to bam files
    # Peaks: path to narrowPeak files (MACS2 output)
    # PeakCaller: type of peak file used (i.e. narrow)

# merged peaks bed file name
MERGED_BED_FILE <- "<NAME_OF_BED_FILE.bed>"


### 1) Load MACS2 peak data ----------------------------------------------------
dba_obj <- dba(sampleSheet = SAMPLE_SHEET)
cat("Finished loading sample sheet and peak data\n")


### 2) Create merged consensus peak BED file -----------------------------------
merged_peaks <- dba.peakset(
  dba_obj,
  bRetrieve = TRUE
)

# Save as BED file (Chrom, Start, End, Score)
write.table(
  merged_peaks[, 1:3],
  file = MERGED_BED_FILE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

cat("Done making merged bed file\n")
