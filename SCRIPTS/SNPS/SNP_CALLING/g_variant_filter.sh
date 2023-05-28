#!/bin/bash -l
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J vcf_filter

module load bioinfo-tools
#module load java/sun_jdk1.8.0_151
module load bcftools
module load python/3.8.7


sample=LsinapisSweM-Swe-sin_final_snps

cp /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/"$sample".vcf $SNIC_TMP

cd $SNIC_TMP


bcftools filter \
	-e 'AF<0.2 || INFO/DP<62 || INFO/DP>194 || QD<2 || FS>60 || MQ<40' \
	-o "$sample"-filtered.vcf \
	"$sample".vcf

cp "$sample"-filtered.vcf /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/


#DP cutoffs: mean +/- 1.5 x sd
#DP mean=128 SD=44 -> cutoffs= 128 + 66 = 194, 128 - 66 = 62


sample=LsinapisSweM-Ire-juv_final_snps

cp /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/"$sample".vcf $SNIC_TMP

cd $SNIC_TMP


bcftools filter \
        -e 'AF<0.2 || INFO/DP<62 || INFO/DP>194 || QD<2 || FS>60 || MQ<40' \
        -o "$sample"-filtered.vcf \
        "$sample".vcf

cp "$sample"-filtered.vcf /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/
