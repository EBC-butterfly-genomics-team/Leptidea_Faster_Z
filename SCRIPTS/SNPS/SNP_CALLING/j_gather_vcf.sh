#!/bin/bash -l
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH -J gather_vcf

module load bioinfo-tools
module load GATK/4.1.1.0
module load java/sun_jdk1.8.0_151


cd /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2

for i in \
	Swe-sin-101C_AGTTCC \
	Swe-sin-102C_ATGTCA \
	Swe-sin-1C_GCCAAT \
	Swe-sin-2C_CAGATC \
	Swe-sin-31C_ACTTGA \
	Swe-sin-32C_GATCAG \
	Swe-sin-61C_TAGCTT \
	Swe-sin-62C_GGCTAC \
	Swe-sin-91C_CTTGTA \
	Swe-sin-92C_AGTCAA \
	Ire-juv-1C_CCGTCC \
	Ire-juv-21C_GTGAAA \
	Ire-juv-22C_GTGGCC \
	Ire-juv-2C_GTCCGC \
	Ire-juv-41C_GTTTCG \
	Ire-juv-42C_CGTACG \
	Ire-juv-61C_GAGTGG \
	Ire-juv-62C_ACTGAT \
	Ire-juv-81C_ATTCCT \
	Ire-juv-82C_ATCACG

do

  chr_vcf=LsinapisSweM-"$i"-sorted.deduped.bqsr

  gatk --java-options "-Xmx6g" GatherVcfs \
	-I "$chr_vcf"_HiC_scaffold_1.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_2.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_3.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_4.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_5.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_6.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_7.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_8.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_9.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_10.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_11.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_12.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_13.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_14.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_15.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_16.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_17.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_18.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_19.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_20.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_21.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_22.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_23.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_24.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_25.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_26.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_27.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_28.g.vcf.gz \
	-I "$chr_vcf"_HiC_scaffold_29.g.vcf.gz \
	-O "$chr_vcf"_all.g.vcf.gz

  gatk --java-options "-Xmx6g" IndexFeatureFile \
        -F "$chr_vcf"_all.g.vcf.gz \

done
