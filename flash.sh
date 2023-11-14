#!/bin/bash

# This script uses flash to combine short overlapping paired reads into unpaired reads. 
# As minimum requirement an overlap of 10 nucleotides is chosen. 
# Note that flash somehow doesn't take absolute paths as output, for that reason path is defined
# relative. 

#SBATCH -N 1 
#SBATCH --mem=90000 
#SBATCH -o /pfs/work7/workspace/scratch/hd_kg304-YZ619/12_Slurm_reports/flash.%A.out
#SBATCH -p single
#SBATCH -t 02:00:00

declare -a e_array=(0 1)
for e in ${e_array[@]}; do 

	input_dir=/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/02_trimmed_fastq/e_$e
	cd $input_dir
	declare -a treatments=("Control" "FK866" "NRH" "Rotenone")

	for i in {03..08}; do
	    for treat in "${treatments[@]}"; do
		echo "paired_trimmed_bc${i}_${treat}_R1.fastq and paired_trimmed_bc${i}_${treat}_R2.fastq"

		flash \
		    -m 10 \
		    -M 150 \
		    -o "../../03_merged_fastq/e_$e/merged_trimmed_bc${i}_${treat}.fastq" \
		    "paired_trimmed_bc${i}_${treat}_R1.fastq" "paired_trimmed_bc${i}_${treat}_R2.fastq"
	    done
	done
done
