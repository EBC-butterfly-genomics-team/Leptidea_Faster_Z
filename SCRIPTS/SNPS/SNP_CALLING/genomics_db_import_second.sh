#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J genomics_db_import


module load bioinfo-tools
module load GATK/4.1.1.0
module load java/sun_jdk1.8.0_151


cd /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2


gatk --java-options "-Xmx6g" GenomicsDBImport \
	-V LsinapisSweM-Swe-sin-101C_AGTTCC-sorted.deduped.bqsr_all.g.vcf.gz \
	-V LsinapisSweM-Swe-sin-31C_ACTTGA-sorted.deduped.bqsr_all.g.vcf.gz \
	-V LsinapisSweM-Swe-sin-91C_CTTGTA-sorted.deduped.bqsr_all.g.vcf.gz \
	-V LsinapisSweM-Swe-sin-102C_ATGTCA-sorted.deduped.bqsr_all.g.vcf.gz \
	-V LsinapisSweM-Swe-sin-32C_GATCAG-sorted.deduped.bqsr_all.g.vcf.gz \
	-V LsinapisSweM-Swe-sin-92C_AGTCAA-sorted.deduped.bqsr_all.g.vcf.gz \
	-V LsinapisSweM-Swe-sin-1C_GCCAAT-sorted.deduped.bqsr_all.g.vcf.gz \
	-V LsinapisSweM-Swe-sin-61C_TAGCTT-sorted.deduped.bqsr_all.g.vcf.gz \
	-V LsinapisSweM-Swe-sin-2C_CAGATC-sorted.deduped.bqsr_all.g.vcf.gz \
	-V LsinapisSweM-Swe-sin-62C_GGCTAC-sorted.deduped.bqsr_all.g.vcf.gz \
	--genomicsdb-workspace-path LsinapisSweM-Swe-sin_HiC_scaffold_"$1" \
	-L HiC_scaffold_"$1"

gatk --java-options "-Xmx6g" GenomicsDBImport \
	-V LsinapisSweM-Swe-sin-101C_AGTTCC-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Swe-sin-31C_ACTTGA-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Swe-sin-91C_CTTGTA-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Swe-sin-102C_ATGTCA-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Swe-sin-32C_GATCAG-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Swe-sin-92C_AGTCAA-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Swe-sin-1C_GCCAAT-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Swe-sin-61C_TAGCTT-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Swe-sin-2C_CAGATC-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Swe-sin-62C_GGCTAC-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Ire-juv-1C_CCGTCC-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Ire-juv-21C_GTGAAA-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Ire-juv-22C_GTGGCC-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Ire-juv-2C_GTCCGC-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Ire-juv-41C_GTTTCG-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Ire-juv-42C_CGTACG-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Ire-juv-61C_GAGTGG-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Ire-juv-62C_ACTGAT-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Ire-juv-81C_ATTCCT-sorted.deduped.bqsr_all.g.vcf.gz \
        -V LsinapisSweM-Ire-juv-82C_ATCACG-sorted.deduped.bqsr_all.g.vcf.gz \
        --genomicsdb-workspace-path LsinapisSweM-sin-juv_HiC_scaffold_"$1" \
        -L HiC_scaffold_"$1"
