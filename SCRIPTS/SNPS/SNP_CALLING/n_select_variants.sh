#!/bin/bash -l
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J select_snps

module load bioinfo-tools
module load GATK/4.1.1.0
module load java/sun_jdk1.8.0_151
module load samtools


reference=LsinapisSweM.fasta
sinapis=LsinapisSweM-Swe-sin
sinjuv=LsinapisSweM-sin-juv


cp /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/$reference $SNIC_TMP


cd $SNIC_TMP

samtools faidx $reference
gatk --java-options "-Xmx6g" CreateSequenceDictionary -R $reference

# sin only...

cp /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2/"$sinapis"_final_all.vcf* $SNIC_TMP

gatk --java-options "-Xmx6g" IndexFeatureFile \
        -F "$sinapis"_final.vcf.gz


gatk --java-options "-Xmx6g" SelectVariants \
	-R $reference \
	-V "$sinapis"_final.vcf.gz \
	--select-type-to-include SNP \
	-O "$sinapis"_final-snps.vcf.gz

cp "$sinapis"_final-snps.vcf.* /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2


# sin+juv...

cp /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2/"$sinjuv"_final_all.vcf* $SNIC_TMP

gatk --java-options "-Xmx6g" IndexFeatureFile \
        -F "$sinjuv"_final_all.vcf.gz


gatk --java-options "-Xmx6g" SelectVariants \
        -R $reference \
        -V "$sinjuv"_final_all.vcf.gz \
        --select-type-to-include SNP \
        -O "$sinjuv"_final_all-snps.vcf.gz

cp "$sinjuv"_final_all-snps.vcf.* /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2
