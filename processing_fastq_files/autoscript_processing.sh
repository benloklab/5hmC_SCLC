#!/bin/bash
#SBATCH -t 00:30:00
#SBATCH -p all
#SBATCH --mem=256M
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o %x-%j.out

##################################################################################################
### Automation script
# This script makes custom scripts for each sample to process fastq files and align them to the human genome
	# Uses trimgalore!, bwa-mem, samtools
# Creator: Janice Li
# Last edit: Mar 18, 2024
##################################################################################################


### Get all unique file names for which automated scripts should be made (e.g. TH4283B_B11)
#all_files=$(ls *.fastq.gz | cut -d '_' -f 1-2| sort -u)
all_files=$(ls *.fq.gz | cut -d '_' -f 1-2| sort -u)

echo "Directory is: $(pwd)"

### Loop through each fastq file in fileNames

for f in $all_files;
do

        # Print which file that is being worked on
	echo "Working with file id: ${f}"
	echo "Making ${f}.sh"

        # Create script file and add shebang line
        touch ${f}.sh
        echo '#!/bin/bash' > ${f}.sh

        # Add the sbatch parameters
        echo '#SBATCH -t 72:00:00' >> ${f}.sh
        echo '#SBATCH -p himem' >> ${f}.sh
        echo '#SBATCH --mem=60G' >> ${f}.sh
        echo '#SBATCH -c 1' >> ${f}.sh
        echo '#SBATCH -N 1' >> ${f}.sh
        echo '#SBATCH -o %x-%j.out' >> ${f}.sh
        echo '' >> ${f}.sh
		
        
	# Import modules
        echo '### Import modules' >> ${f}.sh
        echo 'module load trim_galore/0.5.0' >> ${f}.sh
        echo 'module load cutadapt/2.5' >> ${f}.sh
        echo 'module load bwa/0.7.15' >> ${f}.sh
        echo 'module load igenome-human/hg19' >> ${f}.sh
        echo 'module load samtools/1.10' >> ${f}.sh
        echo '' >> ${f}.sh
        echo '' >> ${f}.sh
		
		
        # Run TrimGalore!
        echo '### Trim out adaptors and perform fastqc' >> ${f}.sh
        echo "trim_galore --paired --fastqc ${f}_R1.fastq.gz ${f}_R2.fastq.gz -o trimmed_sequences/" >> ${f}.sh
        echo '' >> ${f}.sh
		
		
        # Change directory to the trimmed_sequences directory (sam files in here)
        echo '# Change directory to trimmed_sequences' >> ${f}.sh
		echo 'cd trimmed_sequences/' >> ${f}.sh
		
		### Set variables
		echo '' >> ${f}.sh
		echo '# Set variables (for aligning to human genome)' >> ${f}.sh
		echo "sam_file=${f}.sam" >> ${f}.sh
        echo "seq_1=${f}_R1_val_1.fq.gz" >> ${f}.sh
        echo "seq_2=${f}_R2_val_2.fq.gz" >> ${f}.sh
		echo '' >> ${f}.sh
	
        ### Run bwa mem and align to human genome
        echo '### Align sequences to human genome' >> ${f}.sh
        echo 'bwa mem -M $BWAINDEX $seq_1 $seq_2 >> $sam_file' >> ${f}.sh
		echo  "echo 'Done making sam file'" >> ${f}.sh
		echo '' >> ${f}.sh

		echo "bam_file=${f}.bam" >> ${f}.sh
        echo "bam_file_sorted=${f}.sorted.bam" >> ${f}.sh
		echo "bam_file_rmdup=${f}.rmdp.bam" >> ${f}.sh
			
        
        ### Convert SAM file to BAM file
		echo '' >> ${f}.sh
		echo '### Convert sam to bam' >> ${f}.sh	
	
		# Note that we only removed with -f 2 and -F 4 because further filtering will occur later on
		echo 'samtools view -bS -f 2 -F 4 $sam_file | samtools sort -@4 -o $bam_file' >> ${f}.sh
		echo '' >> ${f}.sh
	
		# Deletes original SAM
        echo '# Remove sam file and other bam files' >> ${f}.sh
        echo 'rm $sam_file' >> ${f}.sh
		echo  "echo 'Deleted sam file'" >> ${f}.sh
		
        ### Sort, remove duplicates, and index bam file
		echo '' >> ${f}.sh
		# Sort 
        echo '### Sort, remove duplicates, and index bam file' >> ${f}.sh	
		echo 'samtools sort $bam_file -o $bam_file_sorted' >> ${f}.sh
        # Remove duplicates
        echo 'samtools rmdup $bam_file_sorted $bam_file_rmdup' >> ${f}.sh
        # Indexing BAM file
        echo 'samtools index $bam_file_rmdup' >> ${f}.sh
		echo '' >> ${f}.sh
		
        # Deletes original BAM (contains duplicate reads)
		echo '# Remove other bam files' >> ${f}.sh
		echo 'rm $bam_file' >> ${f}.sh
		echo 'rm $bam_file_sorted' >> ${f}.sh        

		echo '' >> ${f}.sh
		# Computing flagstat and stats
		echo '### Compute flagstat and stat' >> ${f}.sh
        echo 'samtools flagstat $bam_file_rmdup' >> ${f}.sh
        echo 'samtools stats $bam_file_rmdup' >> ${f}.sh
	
		echo '' >> ${f}.sh
		echo "echo 'Done making bam files aligning to the human genome'" >> ${f}.sh

done

echo "Done making scripts"

