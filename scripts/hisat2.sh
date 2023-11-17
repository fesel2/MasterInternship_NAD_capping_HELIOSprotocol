#!/bin/bash

# This script uses hisat2 to align the reads to the reference genome which was downloaded
# from ensemble.
# Note that it requires building indices beforehand.
# It alignes the merged reads (combined by flash) separately than the unmerged (long reads)

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --export=ALL,OMP_NUM_THREADS=8
#SBATCH --mem=5000
#SBATCH --mail-user=kg304@uni-heidelberg.de 
#SBATCH --mail-type=END
#SBATCH -o /pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/12_Slurm_reports/output.%A.out
#SBATCH -t 06:00:00
#SBATCH -p single

declare -a e_array=(0 1)

for e in ${e_array[@]}; do

	export HISAT2_INDEXES=/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/10_alignment_data/HISAT_indexes

	input_dir=/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/03_merged_fastq/e_$e
	output_dir=/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/04_aligned_sequences/e_$e

	cd $input_dir

	declare -a forward_reads=(*notCombined_1.fastq)
	declare -a reverse_reads=(*notCombined_2.fastq)
	declare -a merged_reads=(*extendedFrags.fastq)

	# Loop for paired Reads after flash
	i=0

	for reads in ${forward_reads[@]}; do

		new_label=$(echo $reads | cut -d . -f 1)
		echo "Aligning Sequences of ${forward_reads[i]}, ${reverse_reads[i]}, saving in name $output_dir/$new_label.paired.sam"

		hisat2 \
			-p 4 \
			-x GRCh38.p14.genome \
			-q \
			-1 ${forward_reads[i]} \
			-2 ${reverse_reads[i]} \
			-S $output_dir/$new_label.paired.sam
		i=$((i+1))
	done

	# Loop for merged Reads
	i=0

	for reads in ${forward_reads[@]}; do

		new_label=$(echo $reads | cut -d . -f 1)
		echo "Aligning Sequences of ${merged_reads[i]}, saving in name $output_dir/$new_label.short.sam"

		hisat2 \
			-p 4 \
			-x GRCh38.p14.genome \
			-q \
			-U ${merged_reads[i]} \
			-S $output_dir/$new_label.short.sam
		i=$((i+1))
	done
done
