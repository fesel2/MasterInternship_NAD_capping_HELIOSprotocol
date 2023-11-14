#!/bin/bash

# this script decompresses the originalfastq files into standard .fastq format
# additionally it runs fastQC on all the extracted files
# after that i concatenates the Lanes and saves them in a new directory 
# multiqc report in parent directory combines result of all the fastQC reports
# it requires the Anaconda environment NADcapturedSeq to be activatet to run
# Before running, original data was unzipped and put into the folder unpacked

##SBATCH -N 1
#SBATCH --mem=90000
#SBATCH --mail-user=kg304@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH -o /pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/Slurm_reports/output.%A.out
#SBATCH -t 10:00:00
#SBATCH -p single

start_dir=/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/00_original_data

cd $start_dir
mkdir concatenated_fastq

cd unpacked
declare -a directories=($(ls))

# iterate over every condition, Read-direction and Lane to gather data and run fastqc
for treat in ${directories[@]}; do 

	cd $treat/${treat}_Read1
	gunzip *.gz
	rm *.gz
	for file in *.fastq; do
		fastqc $file
	done
		
	cat *.fastq > ../../../concatenated_fastq/${treat}_R1.fastq

	cd ../..

	cd $treat/${treat}_Read2
	gunzip *.gz
	rm *.gz
	for file in *.fastq; do
		fastqc $file
	done

	cat *.fastq > ../../../concatenated_fastq/${treat}_R2.fastq

	cd ../..
done

# run multiqc for each sample to have a check the lanes in each sample
multiqc . 

# run fastqc and multiqc for the concatenated fastq files
cd ../concatenated_fastq
for file in *.fastq; do
	fastqc $file
done
multiqc . 

