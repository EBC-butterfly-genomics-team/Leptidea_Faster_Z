#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J gather_all.vcf

module load bioinfo-tools
module load GATK/4.1.1.0
module load java/sun_jdk1.8.0_151


cd /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2


gatk --java-options "-Xmx6g" GatherVcfs \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_1_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_2_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_3_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_4_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_5_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_6_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_7_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_8_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_9_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_10_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_11_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_12_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_13_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_14_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_15_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_16_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_17_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_18_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_19_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_20_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_21_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_22_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_23_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_24_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_25_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_26_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_27_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_28_all.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_29_all.vcf \
	 -O LsinapisSweM-Swe-sin_final_all.vcf.gz

gatk --java-options "-Xmx6g" IndexFeatureFile \
        -F LsinapisSweM-Swe-sin_final_all.vcf.gz



gatk --java-options "-Xmx6g" GatherVcfs \
        -I LsinapisSweM-sin-juv_HiC_scaffold_1_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_2_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_3_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_4_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_5_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_6_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_7_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_8_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_9_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_10_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_11_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_12_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_13_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_14_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_15_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_16_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_17_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_18_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_19_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_20_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_21_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_22_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_23_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_24_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_25_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_26_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_27_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_28_all.vcf \
        -I LsinapisSweM-sin-juv_HiC_scaffold_29_all.vcf \
        -O LsinapisSweM-sin-juv_final_all.vcf.gz

gatk --java-options "-Xmx6g" IndexFeatureFile \
        -F LsinapisSweM-sin-juv_final_all.vcf.gz
