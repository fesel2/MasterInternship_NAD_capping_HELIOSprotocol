#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH --mail-user=kg304@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH -o ./Slurm_reports/output.%A.out
#SBATCH -t 0:10:00
#SBATCH -p single

dir=/pfs/work7/workspace/scratch/hd_kg304-YZ619/01_trim3_demultiplexed_fastq
declare -a e_array=(e_0.03 e_0.06 e_0.09)
for e in ${e_array[@]}; do
	cd $dir/$e 
	declare -a conditions=(Control FK866 Rotenone NRH)

	for cond in ${conditions[@]}; do
		
		reads=$(echo "$(cat *$cond*1*.fastq | wc -l) /4" | bc)
		echo "$cond has $reads reads." >> readcount.txt
	done
done
