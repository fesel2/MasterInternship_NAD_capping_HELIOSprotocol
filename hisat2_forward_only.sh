#!/bin/bash

# This script aims to answer the question if alignment rates of the spikRNA are changing using
# only the forward read.

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --export=ALL,OMP_NUM_THREADS=8
#SBATCH --mem=5000
#SBATCH --mail-user=kg304@uni-heidelberg.de #SBATCH --mail-type=END
#SBATCH -o /pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/12_Slurm_reports/forward_align.%A.out
#SBATCH -t 06:00:00
#SBATCH -p single

declare -a e_array=(1)

for e in ${e_array[@]}; do

	export HISAT2_INDEXES=/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/10_alignment_data/HISAT_indexes

	input_dir=/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/03_merged_fastq/e_$e
	output_dir=/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/04_aligned_sequences/forward_align

	cd $input_dir

	declare -a forward_reads=(*notCombined_1.fastq)
	declare -a reverse_reads=(*notCombined_2.fastq)
	declare -a merged_reads=(*extendedFrags.fastq)

	# Loop forward readsReads

	for reads in ${forward_reads[@]}; do

		new_label=$(echo $reads | cut -d . -f 1)
		echo "Aligning Sequences of $reads (only forward), saving in name $output_dir/$new_label.short.sam"

		hisat2 \
			-p 4 \
			-x GRCh38.p14.genome \
			-q \
			-U $reads \
			-S $output_dir/$new_label.forward.sam
	done

done
