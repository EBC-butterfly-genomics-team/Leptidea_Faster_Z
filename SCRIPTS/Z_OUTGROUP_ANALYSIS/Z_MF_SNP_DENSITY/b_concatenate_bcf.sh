#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J concatenate_bcf

module load bioinfo-tools
module load bcftools

cd /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/SNPS/MPILEUP

for sample in P14458_103_S24 P14458_104_S25 P14458_107_S28 P14458_108_S29 P14502_103 P14502_104
do
  bcftools concat \
	-Ob \
	-o "$sample".bcf \
	"$sample"-HiC_scaffold_1.bcf \
	"$sample"-HiC_scaffold_2.bcf \
	"$sample"-HiC_scaffold_3.bcf \
	"$sample"-HiC_scaffold_4.bcf \
	"$sample"-HiC_scaffold_5.bcf \
	"$sample"-HiC_scaffold_6.bcf \
	"$sample"-HiC_scaffold_7.bcf \
	"$sample"-HiC_scaffold_8.bcf \
	"$sample"-HiC_scaffold_9.bcf \
	"$sample"-HiC_scaffold_10.bcf \
	"$sample"-HiC_scaffold_11.bcf \
	"$sample"-HiC_scaffold_12.bcf \
	"$sample"-HiC_scaffold_13.bcf \
	"$sample"-HiC_scaffold_14.bcf \
	"$sample"-HiC_scaffold_15.bcf \
	"$sample"-HiC_scaffold_16.bcf \
	"$sample"-HiC_scaffold_17.bcf \
	"$sample"-HiC_scaffold_18.bcf \
	"$sample"-HiC_scaffold_19.bcf \
	"$sample"-HiC_scaffold_20.bcf \
	"$sample"-HiC_scaffold_21.bcf \
	"$sample"-HiC_scaffold_22.bcf \
	"$sample"-HiC_scaffold_23.bcf \
	"$sample"-HiC_scaffold_24.bcf \
	"$sample"-HiC_scaffold_25.bcf \
	"$sample"-HiC_scaffold_26.bcf \
	"$sample"-HiC_scaffold_27.bcf \
	"$sample"-HiC_scaffold_28.bcf \
	"$sample"-HiC_scaffold_29.bcf
done
