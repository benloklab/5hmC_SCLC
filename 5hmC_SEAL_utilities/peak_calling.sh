#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -p himem
#SBATCH --mem=60G
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o %x-%j.out

##############################################################################
#
# This script was used to profile 5hmC peaks in an HPC environment
# 
# MACS2
# `INPUT: bam files 
# `OUTPUT: narrowpeak files (OPTIONAL: can choose to output as bed files)
#
# Post-processing:
#   1. Remove ENCODE blacklist regions
#   2. Remove chrX and chrY peaks
#
# Author: Janice Li
# Last updated: May 21, 2026
#
##############################################################################

set -euo pipefail

### Load modules -------------------------------------------------------------
module load gcc/6.2.0
module load python/3.4.3
#module load MACS/2.1.2
module load R/4.2.1
module load bedtools

### Set up -------------------------------------------------------------
BAM_FILE_DIR="<FOLDER_WITH_BAM_FILES>"
MACS2_OUTPUT_DIR="<OUTPUT_FOLDER>"
BLACKLIST="<PATH_TO_ENCODE_BLACKLIST_BED_FILE>"
FILTERED_DIR="${MACS2_OUTPUT_DIR}/filtered_narrowPeak"

# Make folders if not already made
#mkdir -p "$MACS2_OUTPUT_DIR"
#mkdir -p "$FILTERED_DIR"

### 1) Perform MACS2 on BAMS ------------------------------------------------------------------------
echo "MACS2 started at: $(date)"
echo "Starting MACS2 peak calling"

# Get all bam files
for BAM_FILE in "${BAM_FILE_DIR}"/*.bam; do
  SAMPLE_NAME="$(basename "$BAM_FILE" .bam)"
  SAMPLE_NAME="${SAMPLE_NAME%.rmdp}"
  
  echo "Running MACS2 for: $BAM_FILE"
  
  # Run MACS2
   macs2 callpeak -t "$BAM_FILE" -f BAM -p 1e-5 -n "$SAMPLE_NAME" --outdir "$MACS2_OUTPUT_DIR"

  echo "Finished peak calling for: $SAMPLE_NAME"
done

echo "MACS2 finished at: $(date)"

### 2) Filter narrowPeak files ------------------------------------------------------------------------
# Remove ENCODE blacklist regions and X/Y chromosomes 
# Note: if bed file outputs are desired, replace "narrowPeak" with "bed"

echo "Filtering started at: $(date)"
echo "Starting narrowPeak filtering"

for PEAK_FILE in "${MACS2_OUTPUT_DIR}"/*.narrowPeak; do
  # Define output file names
  PEAK_BASENAME="$(basename "$PEAK_FILE")"
  SAMPLE_NAME="${PEAK_BASENAME%_peaks.narrowPeak}"

  FILTERED_PEAKS="${FILTERED_DIR}/${SAMPLE_NAME}.bf.narrowPeak"
  FINAL_PEAKS="${FILTERED_DIR}/${SAMPLE_NAME}.bf_noXY.narrowPeak"

  echo "Filtering: $PEAK_BASENAME"
  
  # Remove ENCODE blacklist regions
  bedtools intersect -v -a "$PEAK_FILE" -b "$BLACKLIST" > "$FILTERED_PEAKS"
  echo "Removed blacklist regions from $PEAK_FILE"
  
  # Remove sex chromosomes
  awk '$1 != "chrX" && $1 != "chrY"' "$FILTERED_PEAKS" > "$FINAL_PEAKS"
  echo "Removed X/Y chromosomes from $PEAK_FILE"
 
done

echo "Filtering narrowPeaks files complete"
echo "Filtering finished at: $(date)"
