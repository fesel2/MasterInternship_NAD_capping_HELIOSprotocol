#!/bin/bash

# This script alignes the reads to the known sequences of the spikeRNAs. 
# Make sure to have built the hisat indices beforehand. Each of the indices requires a separate 
# folder to work. then this script can loop over the indeces

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --export=ALL,OMP_NUM_THREADS=8
#SBATCH --mem=5000
#SBATCH --mail-user=kg304@uni-heidelberg.de #SBATCH --mail-type=END
#SBATCH -o /pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/12_Slurm_reports/hisat2_spikeRNA_align.%A.out
#SBATCH -t 06:00:00
#SBATCH -p single


index_dir=/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/10_alignment_data/HISAT_indexes
start_dir=/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/03_merged_fastq/e_1
output_dir=/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/04_aligned_sequences/spike_alignment

declare -a spike_list=(biotinRna nad36 nad60 nad102A nad200 nad351 nad501 m7G122 m7G250 m7G400 ppp301 ppp51 ppp102 ppp149)
cd $start_dir

declare -a forward_reads=(*notCombined_1.fastq)
declare -a reverse_reads=(*notCombined_2.fastq)
declare -a merged_reads=(*extendedFrags.fastq)


for spike in ${spike_list[@]}; do 
	export HISAT2_INDEXES=/pfs/work7/workspace/scratch/hd_kg304-YZ619/10_alignment_data/HISAT_indexes/$spike

	# Loop for paired Reads after flash
        i=0

        for reads in ${forward_reads[@]}; do

                new_label=$(echo $reads | cut -d . -f 1)
                echo "Aligning Sequences of ${forward_reads[i]}, ${reverse_reads[i]} to $spike sequence"

                hisat2 \
                        -p 4 \
                        -x $spike \
                        -q \
                        -1 ${forward_reads[i]} \
                        -2 ${reverse_reads[i]} \
                        -S $output_dir/$spike/$new_label.paired.sam
                i=$((i+1))
	done

	# Loop for merged Reads
        i=0

        for reads in ${forward_reads[@]}; do

                new_label=$(echo $reads | cut -d . -f 1)
                echo "Aligning Sequences of ${merged_reads[i]} to $spike sequence."

                hisat2 \
                        -p 4 \
                        -x $spike \
                        -q \
                        -U ${merged_reads[i]} \
                        -S $output_dir/$spike/$new_label.short.sam
                i=$((i+1))
        done
done
