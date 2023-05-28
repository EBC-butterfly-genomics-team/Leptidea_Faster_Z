#!/bin/bash -l
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J gather_vcf

module load bioinfo-tools
module load GATK/4.1.1.0
module load java/sun_jdk1.8.0_151


cd /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller

gatk --java-options "-Xmx6g" GatherVcfs \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_1.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_2.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_3.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_4.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_5.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_6.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_7.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_8.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_9.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_10.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_11.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_12.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_13.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_14.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_15.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_16.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_17.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_18.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_19.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_20.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_21.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_22.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_23.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_24.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_25.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_26.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_27.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_28.vcf \
        -I LsinapisSweM-Swe-sin_HiC_scaffold_29.vcf \
	-O LsinapisSweM-Swe-sin_final.vcf.gz

gatk --java-options "-Xmx6g" IndexFeatureFile \
        -F LsinapisSweM-Swe-sin_final.vcf.gz

gatk --java-options "-Xmx6g" GatherVcfs \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_1.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_2.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_3.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_4.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_5.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_6.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_7.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_8.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_9.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_10.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_11.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_12.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_13.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_14.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_15.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_16.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_17.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_18.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_19.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_20.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_21.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_22.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_23.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_24.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_25.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_26.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_27.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_28.vcf \
        -I LsinapisSweM-Ire-juv_HiC_scaffold_29.vcf \
	-O LsinapisSweM-Ire-juv_final.vcf.gz

gatk --java-options "-Xmx6g" IndexFeatureFile \
        -F LsinapisSweM-Ire-juv_final.vcf.gz
