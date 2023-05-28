#!/bin/bash -l
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J select_snps

module load bioinfo-tools
module load GATK/4.1.1.0
module load java/sun_jdk1.8.0_151
module load samtools


reference=LsinapisSweM.fasta


for sample in LsinapisSweM-Swe-sin LsinapisSweM-Ire-juv
do

  cp /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/$reference $SNIC_TMP
  cp /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/"$sample"_final.vcf.* $SNIC_TMP

  cd $SNIC_TMP

  samtools faidx $reference
  gatk --java-options "-Xmx6g" CreateSequenceDictionary -R $reference

  gatk --java-options "-Xmx6g" SelectVariants \
	-R $reference \
	-V "$sample"_final.vcf.gz \
	--select-type-to-include SNP \
	-O "$sample"_final_snps.vcf

  cp "$sample"_final_snps.vcf /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller

done
