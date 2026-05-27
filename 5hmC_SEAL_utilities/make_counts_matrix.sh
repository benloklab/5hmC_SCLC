#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -p veryhimem
#SBATCH --mem=400G
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -o %x-%j.out

##############################################################################
#
# This script was used to run DiffBind in R and featureCounts in an HPC environment
# 
# The setup above can be modified but note that Diffbind is a memory-heavy
# process, especially with >100 samples. If just running featureCounts, 32G 
# for 12:00:00 may suffice
# 
# Diffbind
# 	INPUT: csv file with column headers: Sample_ID, Condition, Replicate, bamReads, Peak, and Peakcaller
# 	OUTPUT: BED file 
#
# featureCounts
# 	INPUTS:
#		cons: bed file containing merged peakset (output from diffbind.R)
#		bams: text file containing the paths to all bam files of interest
# 	OUTPUT: 
#		counts.txt: counts data
#		counts_mat.txt: counts matrix 
#	
#
# Author Janice Li
# Last updated: May 21, 2026
#
##############################################################################

set -euo pipefail

### Load modules -------------------------------------------------------------
module load R/4.2.1
module load subread

# User-defined paths -----------------------------------------------------------
CONSENSUS_PEAKS="merged_peaks.bed"
BAM_LIST="bam_files.txt"
OUTPUT_DIR="<OUTPUT_DIRECTORY_NAME>"

mkdir -p "$OUTPUT_DIR"

### Run diffbind.R to make merged peak set ----------------------------------------
echo 'Starting R driver script (R version 4.2.1)'
echo "Diffbind started at:" $(date)

Rscript diffbind.R

echo 'Finished running diffbind.R'
echo "Diffbind finished at:" $(date)


### Run featureCounts ------------------------------------------------------------

echo "FeatureCounts started at: $(date)"
# Prepare SAF annotation file
echo "Making saf file"

SAF_FILE="${CONSENSUS_PEAKS%.bed}.saf"

echo "Creating SAF file: $SAF_FILE"

awk 'BEGIN {OFS="\t"} {print $1"."$2"."$3, $1, $2, $3, "."}' \
  "$CONSENSUS_PEAKS" > "$SAF_FILE"

# Run featureCounts
echo "Running feature counts"
featureCounts -a $SAF_FILE -o $OUTPUT_DIR/counts.txt -F 'SAF' --primary -T 8 -p $(cat $bams | xargs) 
# Since we work with the filtered files, we do not need to worry about duplicates or bam quality reads

# Extract counts matrix
echo "Extracting counts matrix"
cut -f1,7- $OUTPUT_DIR/counts.txt > $OUTPUT_DIR/counts_mat.txt

echo "Done making counts matrix"
echo "FeatureCounts finished at:" $(date)
