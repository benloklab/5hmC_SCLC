#!/bin/bash
#SBATCH -t 00:30:00
#SBATCH -p all
#SBATCH --mem=256M
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o %x-%j.out

##################################################################################################
### Automation script
# This script makes custom scripts to map to spike-in controls for each sample
# Runs trimgalore!, bwa-mem, samtools
# Creator: Janice Li
# Last edit: Jun 18, 2024
##################################################################################################


### Get all unique file names for which automated scripts should be made
# Outputs TH####X_B##
all_files=$(ls *.fastq.gz | cut -d '_' -f 1-2| sort -u)

echo "Directory is: $(pwd)"

### Loop through each fastq file

for f in $all_files;
do

	# Print which file that is being worked on
	echo "Working with file id: ${f}"
	echo "Making ${f}.sh"

	# Create script file and add shebang line
	touch ${f}_spike.sh
	echo '#!/bin/bash' > ${f}_spike.sh

	# Add the sbatch parameters
	echo '#SBATCH -t 24:00:00' >> ${f}_spike.sh
	echo '#SBATCH -p himem' >> ${f}_spike.sh
	echo '#SBATCH --mem=60G' >> ${f}_spike.sh
	echo '#SBATCH -c 1' >> ${f}_spike.sh
	echo '#SBATCH -N 1' >> ${f}_spike.sh
	echo '#SBATCH -o %x-%j.out' >> ${f}_spike.sh

	echo '' >> ${f}_spike.sh
		
        
	# Import modules
	echo '### Import modules' >> ${f}_spike.sh
	echo 'module load trim_galore' >> ${f}_spike.sh
	echo 'module load cutadapt' >> ${f}_spike.sh
	echo 'module load bwa' >> ${f}_spike.sh
	echo 'module load igenome-human/hg19' >> ${f}_spike.sh
	echo 'module load samtools' >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh
	
	
	# Run TrimGalore!
	echo '### Trim out adaptors + fastqc' >> ${f}_spike.sh
	echo "# Note that this should already have been done in the preprocessing step so commented out"
	echo "#trim_galore --paired --fastqc ${f}_R1.fastq.gz ${f}_R2.fastq.gz -o <FOLDER_WITH_FASTQ_FILES>" >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh
	
	# Change directory to the trimmed_sequences directory (sam files in here)
	echo '# Change directory to trimmed_spikes' >> ${f}_spike.sh
	echo 'cd trimmed_sequences/' >> ${f}_spike.sh
	
	### Set variables
	echo '' >> ${f}_spike.sh
	echo '# Set variables (for aligning to escherichia phage lambda genome)' >> ${f}_spike.sh
	echo "sam_spike=${f}_spike.sam" >> ${f}_spike.sh
	echo "sam_c=${f}_c.sam" >> ${f}_spike.sh
	echo "sam_m=${f}_m.sam" >> ${f}_spike.sh
	echo "sam_h=${f}_h.sam" >> ${f}_spike.sh 
	echo "seq_1=${f}_R1_val_1.fq.gz" >> ${f}_spike.sh
	echo "seq_2=${f}_R2_val_2.fq.gz" >> ${f}_spike.sh
	echo "spike_ref=<FASTA_REF_FILE_LOCATION>" >> ${f}_spike.sh
	echo "c_spike_ref=<C_SPIKE_FASTA_REF_FILE_LOCATION>" >> ${f}_spike.sh
	echo "m_spike_ref=<M_SPIKE_FASTA_REF_FILE_LOCATION>" >> ${f}_spike.sh
	echo "h_spike_ref=<H_SPIKE_FASTA_REF_FILE_LOCATION>" >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh

	### Run bwa mem and align to human genome
	echo '### Align sequences to human genome' >> ${f}_spike.sh
	echo 'bwa mem -M $spike_ref $seq_1 $seq_2 >> $sam_spike' >> ${f}_spike.sh
	echo 'bwa mem -M $c_spike_ref $seq_1 $seq_2 >> $sam_c' >> ${f}_spike.sh
	echo 'bwa mem -M $m_spike_ref $seq_1 $seq_2 >> $sam_m' >> ${f}_spike.sh
	echo 'bwa mem -M $h_spike_ref $seq_1 $seq_2 >> $sam_h' >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh
	
	# bam variables for spike (all)
	echo '# BAM variables for spike (all)' >> ${f}_spike.sh
	echo "bam_file_spike=${f}_spike.bam" >> ${f}_spike.sh
    echo "bam_file_sorted_spike=${f}_spike.sorted.bam" >> ${f}_spike.sh
    echo "bam_file_rmdup_spike=${f}_spike.rmdp.bam" >> ${f}_spike.sh
    echo '' >> ${f}_spike.sh

	# bam variables for C spike
	echo '# BAM variables for C Spike' >> ${f}_spike.sh
	echo "bam_file_c=${f}_c.bam" >> ${f}_spike.sh
	echo "bam_file_sorted_c=${f}_c.sorted.bam" >> ${f}_spike.sh
	echo "bam_file_rmdup_c=${f}_c.rmdp.bam" >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh
	
	# bam variables for M spike
	echo '# BAM variables for M Spike' >> ${f}_spike.sh
	echo "bam_file_m=${f}_m.bam" >> ${f}_spike.sh
	echo "bam_file_sorted_m=${f}_m.sorted.bam" >> ${f}_spike.sh
	echo "bam_file_rmdup_m=${f}_m.rmdp.bam" >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh
	
	# bam variables for H spike
	echo '# BAM variables for H Spike' >> ${f}_spike.sh
	echo "bam_file_h=${f}_h.bam" >> ${f}_spike.sh
	echo "bam_file_sorted_h=${f}_h.sorted.bam" >> ${f}_spike.sh
	echo "bam_file_rmdup_h=${f}_h.rmdp.bam" >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh
	
	### Convert SAM file to BAM file
	echo '' >> ${f}_spike.sh
	echo '### Convert sam to bam' >> ${f}_spike.sh

	# Note that we only removed with -f 2 and -F 4 because further filtering will occur later on
	echo 'samtools view -bS -f 2 -F 4 $sam_spike | samtools sort -@4 -o $bam_file_spike' >> ${f}_spike.sh
	echo 'samtools view -bS -f 2 -F 4 $sam_c | samtools sort -@4 -o $bam_file_c' >> ${f}_spike.sh
	echo 'samtools view -bS -f 2 -F 4 $sam_m | samtools sort -@4 -o $bam_file_m' >> ${f}_spike.sh
	echo 'samtools view -bS -f 2 -F 4 $sam_h | samtools sort -@4 -o $bam_file_h' >> ${f}_spike.sh

	# Remove SAM
	echo '# Remove SAM' >> ${f}_spike.sh
	echo 'rm $sam_spike' >> ${f}_spike.sh
	echo 'rm $sam_c' >> ${f}_spike.sh
    echo 'rm $sam_m' >> ${f}_spike.sh
    echo 'rm $sam_h' >> ${f}_spike.sh
	
	### Sort, remove duplicates, and index bam file
	echo '' >> ${f}_spike.sh
	
	# Sort
	echo '### Sort, remove duplicates, and index bam file' >> ${f}_spike.sh
	echo 'samtools sort $bam_file_spike -o $bam_file_sorted_spike' >> ${f}_spike.sh
	echo 'samtools sort $bam_file_c -o $bam_file_sorted_c' >> ${f}_spike.sh
	echo 'samtools sort $bam_file_m -o $bam_file_sorted_m' >> ${f}_spike.sh
	echo 'samtools sort $bam_file_h -o $bam_file_sorted_h' >> ${f}_spike.sh
	
	# Remove duplicates
	echo 'samtools rmdup $bam_file_sorted_spike $bam_file_rmdup_spike' >> ${f}_spike.sh
	echo 'samtools rmdup $bam_file_sorted_c $bam_file_rmdup_c' >> ${f}_spike.sh
	echo 'samtools rmdup $bam_file_sorted_m $bam_file_rmdup_m' >> ${f}_spike.sh
	echo 'samtools rmdup $bam_file_sorted_h $bam_file_rmdup_h' >> ${f}_spike.sh
	
	# Indexing BAM file
	echo 'samtools index $bam_file_rmdup_spike' >> ${f}_spike.sh
	echo 'samtools index $bam_file_rmdup_c' >> ${f}_spike.sh
	echo 'samtools index $bam_file_rmdup_m' >> ${f}_spike.sh
	echo 'samtools index $bam_file_rmdup_h' >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh

	# Deletes original BAM and SAM (contains duplicate reads)
	echo '# Remove other bam files' >> ${f}_spike.sh
	echo 'rm $bam_file_spike' >> ${f}_spike.sh
	echo 'rm $bam_file_c' >> ${f}_spike.sh
	echo 'rm $bam_file_m' >> ${f}_spike.sh
	echo 'rm $bam_file_h' >> ${f}_spike.sh
	echo 'rm $bam_file_sorted_spike' >> ${f}_spike.sh
	echo 'rm $bam_file_sorted_c' >> ${f}_spike.sh
	echo 'rm $bam_file_sorted_m' >> ${f}_spike.sh
	echo 'rm $bam_file_sorted_h' >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh


	# Computing flagstat and stats
	echo '### Compute flagestat and stat' >> ${f}_spike.sh
	echo 'samtools flagstat $bam_file_rmdup_spike' >> ${f}_spike.sh
    echo 'samtools stats $bam_file_rmdup_spike' >> ${f}_spike.sh
    echo '' >> ${f}_spike.sh
	
	# Computing flagstat and stats
	echo 'samtools flagstat $bam_file_rmdup_c' >> ${f}_spike.sh
	echo 'samtools stats $bam_file_rmdup_c' >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh
	
	echo 'samtools flagstat $bam_file_rmdup_m' >> ${f}_spike.sh
	echo 'samtools stats $bam_file_rmdup_m' >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh
	
	echo 'samtools flagstat $bam_file_rmdup_h' >> ${f}_spike.sh
	echo 'samtools stats $bam_file_rmdup_h' >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh
	echo '' >> ${f}_spike.sh
	
	echo "echo 'Done making bam files aligning to the phage genome'" >> ${f}_spike.sh

done

echo "Done making spike in scripts"

