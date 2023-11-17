#!/bin/bash

# This script uses the concatenated .fastq as input. 
# It demultiplexed based on the 5' end of the reverse read.
# Note that for successfull demultiplexing the original adapter sequence is reverse complementary.
# On the other side it also trims the adapter sequence at the 3' end of the forward read
# using wildcards for the different barcodes.
# Additionally it repeats the process for e=0 and e=1. e is the allowed number of missmatches for 
# demultiplexing. Lower e results in less demultiplexed reads and more reads in unknown category.
# After cutadapt FastQC and multiqc is run to check for adapter content.


#SBATCH -N 1 
#SBATCH --mail-user=kg304@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH -o /pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/12_Slurm_reports/cutadapt_demultiplexing.%A.out
#SBATCH -t 06:00:00
#SBATCH -p single
#SBATCH -c 40 
#SBATCH --export=ALL,OMP_NUM_THREADS=20

start_dir=/pfs/work7/workspace/scratch/hd_kg304-YZ619/00_original_data/concatenated_fastq 
end_dir=/pfs/work7/workspace/scratch/hd_kg304-YZ619/01_trim5_demultiplexed_fastq
cd $start_dir

forward_reads=(*1.fastq)
reverse_reads=(*2.fastq)

i=0

# Analysis is run for 0 and 1 Missmatches for adapters because they all differ in at least 2 bases. 
# 1 might improve Demultiplexing

declare -a e_array=(0 1)

for forward_read in ${forward_reads[@]}; do
        declare -a paired_read=(${forward_reads[$i]} ${reverse_reads[$i]})
	for e in ${e_array[@]}; do

		echo "Starting cutadapt in paired mode on these two files: "${paired_read[0]}" "${paired_read[1]}" \
			with absolute error rate of $e"

		# Note that R1 and R2 need to be switched in order to make demultiplexing with R2 possible
		# The 5' end of the reverse strand is used to demultiplex.
		# Also barcoded adapters are removed at 3' of the forward strand (reverse complements)

		cutadapt \
		-j 0 \
		-g bc03=^ACCGGTNNNNNNG \
		-g bc04=^ATGAGTNNNNNNG \
		-g bc05=^GTTCGTNNNNNNG \
		-g bc06=^TGCTGTNNNNNNG \
		-g bc07=^TATGGTNNNNNNG \
		-g bc08=^CTATGTNNNNNNG \
                -A "CNNNNNNACNNNNAGATCGGAAGAGCACACGTCTG;min_overlap=9" \
		-e $e \
		-o "$end_dir/e_${e}/trimmed_{name}_${paired_read[1]}" \
		-p "$end_dir/e_${e}/trimmed_{name}_${paired_read[0]}" \
		"${paired_read[1]}" "${paired_read[0]}" 
	done

        i=$((i+1))
done

# switch in the subfolders and run quality control
for e in ${e_array[@]};do
	cd "$end_dir/e_${e}"
	for file in *.fastq
		do fastqc $file
	done
	multiqc . 
done
