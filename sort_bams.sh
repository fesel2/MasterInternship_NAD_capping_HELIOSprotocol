#!/bin/bash

# This script uses samtools to sort, index and convert the .sam files for further processing

#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH --mail-user=kg304@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH -o /pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/12_Slurm_reports/output.%A.out
#SBATCH -t 02:00:00
#SBATCH -p single


declare -a e_array=(0 1)
for e in ${e_array[@]};do

	cd /pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/04_aligned_sequences/e_$e

	declare -a filenames=($(ls *.sam | cut -f1-2 -d .))
	i=0
	for file in *.sam; do
		echo "working on $file"
		samtools sort $file -o ${filenames[i]}.sorted.bam
		samtools index $file
		i=$((i+1))

	done
done


