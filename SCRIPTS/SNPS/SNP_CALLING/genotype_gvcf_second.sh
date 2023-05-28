#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J genotype_gvcf


module load bioinfo-tools
module load GATK/4.1.1.0
module load java/sun_jdk1.8.0_151
module load samtools


reference=LsinapisSweM
vcf_path=/home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2

cp /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/$reference.fasta $SNIC_TMP

cp -r $vcf_path/$reference-Swe-sin_HiC_scaffold_"$1" $SNIC_TMP
cp -r $vcf_path/$reference-sin-juv_HiC_scaffold_"$1" $SNIC_TMP

cd $SNIC_TMP

samtools faidx $reference.fasta
gatk --java-options "-Xmx6g" CreateSequenceDictionary -R $reference.fasta



# sinapis

gatk --java-options "-Xmx6g" GenotypeGVCFs \
	-R $reference.fasta \
	-V gendb://$reference-Swe-sin_HiC_scaffold_"$1" \
	-L HiC_scaffold_"$1" \
	-O $reference-Swe-sin_HiC_scaffold_"$1"_all.vcf \
	-all-sites

cp $reference-Swe-sin_HiC_scaffold_"$1"* $vcf_path



# sinapis + juvernica

gatk --java-options "-Xmx6g" GenotypeGVCFs \
	-R $reference.fasta \
	-V gendb://$reference-sin-juv_HiC_scaffold_"$1" \
	-L HiC_scaffold_"$1" \
	-O $reference-sin-juv_HiC_scaffold_"$1"_all.vcf \
	-all-sites

cp $reference-sin-juv_HiC_scaffold_"$1"* $vcf_path
