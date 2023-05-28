#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J vcf_filter

module load bioinfo-tools
module load bcftools
module load python/3.8.7

cd $SNIC_TMP

for sample in LsinapisSweM-Swe-sin_final_all LsinapisSweM-sin-juv_final_all
do

  cp /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2/"$sample".vcf.gz $SNIC_TMP

  # for pixy, set missing as "."...

  bcftools filter \
       -i 'FMT/DP>5 & FMT/DP<25 & QUAL>30' \
       --set-GTs . \
       -o "$sample"-filtered_q30.vcf.gz \
       "$sample".vcf.gz

  cp "$sample"-filtered_q30.vcf.gz /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2

done
