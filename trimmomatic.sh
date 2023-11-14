#!/bin/bash

# This script uses trimmomatic to trim the 5' adapter of the forward read 
# and the 3' adaptor of the reverse read.
# To make that work it requires the sequences supplied in a .fasta file, here in the folder 13_useful stuff
# It also filters low quality reads with the sliding window function 
# and filters reads shorter than 24 nucleotides.

#SBATCH -N 1 
#SBATCH --mail-user=kg304@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH -o /pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/12_Slurm_reports/trimmomatic.%A.out
#SBATCH -t 06:00:00
#SBATCH -p single
#SBATCH -c 40 
#SBATCH --export=ALL,OMP_NUM_THREADS=20

declare -a e_array=(0 1)
for e in ${e_array[@]}; do 
	start_dir=/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/01_trim5_demultiplexed_fastq/e_$e
	end_dir=/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/02_trimmed_fastq/e_$e

	cd $start_dir
	forward_reads=(*1.fastq)
	reverse_reads=(*2.fastq)

	i=0
	for forward_read in ${forward_reads[@]}; do
		declare -a paired_read=(${forward_reads[$i]} ${reverse_reads[$i]})

		echo "starting trimmomatic on ${paired_read[0]} and ${paired_read[1]}"
		trimmomatic \
			PE \
			-threads 20 \
			-phred33 \
			${paired_read[0]} ${paired_read[1]} \
			$end_dir/paired_${paired_read[0]} $end_dir/unpaired_${paired_read[0]} \
			$end_dir/paired_${paired_read[1]} $end_dir/unpaired_${paired_read[1]} \
			ILLUMINACLIP:/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/13_useful_stuff/universal_adaptor.fa:2:30:7 \
			SLIDINGWINDOW:5:20 \
			MINLEN:24
			#LEADING:15 \
			#TRAILING:15 \
			
		i=$((i+1))

	done

	cd $end_dir
	for file in *.fastq; do
		fastqc $file
	done
	multiqc . 
done
