#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J snp_vcf

ml bioinfo-tools
ml bcftools

cd /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2

#sample=LsinapisSweM-Swe-sin_final_all-filtered_q30

sample=LsinapisSweM-Swe-sin_final_all

bcftools view -v snps "$sample".vcf.gz > "$sample"_snps.vcf
