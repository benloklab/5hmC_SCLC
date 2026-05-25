#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -p all
#SBATCH --mem=8G
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o %x-%j.out

##################################################################################################
### Script for combining fastq files from sequencing
# Combine all R1 from the same patient together, combine all R2 from the same patient together
# NOTE: At Sickkids, R1=R1, R2=UMI, R3=R2
# Creator: Janice Li
# Made on: Jun 3, 2024
# Last edit: Jun 6, 2024
##################################################################################################

# Get all fastq files 
# Format: TH + #### + B/R/P + _B## (TH4967B_B08)
all_files=$(ls *.fastq.gz | cut -d '_' -f 1-2| sort -u)

# Combine fastq files of the same T-ID based on R1, R2, or R3
for f in $all_files;
do
	# Combine all R1 together -> read1
	cat $(ls ${f}*R1*.fastq.gz) > ${f}_R1.fastq.gz
	mv ${f}_R1.fastq.gz merged

	# Combine all R2 together -> UMI
	cat $(ls ${f}*R2*.fastq.gz) > ${f}_UMI.fastq.gz
	mv ${f}_UMI.fastq.gz merged

	# Combine all R3 together -> read2
	cat $(ls ${f}*R3*.fastq.gz) > ${f}_R2.fastq.gz
	mv ${f}_R2.fastq.gz merged
done
