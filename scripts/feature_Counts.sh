#!/bin/bash

# This script takes the sorted .bam files output the .gtf annotation file and creates a feature Matrix.
# GTF annotation is from ensemble too. This step must be done twice to avoid confusion.
# First time with the merged (shorter) reads. Second time with the paired (longer) reads. 

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --export=ALL,OMP_NUM_THREADS=8
#SBATCH --mem=50000
#SBATCH --mail-user=kg304@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH -o /pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/12_Slurm_reports/output.%A.out
#SBATCH -t 04:00:00
#SBATCH -p single


declare -a e_array=(0 1)
for e in ${e_array[@]}; do 

	featureCounts -T 8 \
		-p \
		-a /pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/10_alignment_data/gencode.v44.chr_patch_hapl_scaff.annotation.gtf \
		-o /pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/05_count_matrices/e_$e/feature_counts_paired.txt \
		/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/04_aligned_sequences/e_$e/*.paired.sorted.bam

	featureCounts -T 8 \
		-a /pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/10_alignment_data/gencode.v44.chr_patch_hapl_scaff.annotation.gtf \
		-o /pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/05_count_matrices/e_$e/feature_counts_single.txt \
		/pfs/work7/workspace/scratch/hd_kg304-Helios_Hekcells_NADSeq/04_aligned_sequences/e_$e/*.short.sorted.bam
done


