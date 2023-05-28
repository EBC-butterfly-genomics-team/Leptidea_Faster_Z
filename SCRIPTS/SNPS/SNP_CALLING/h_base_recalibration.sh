#!/bin/bash -l
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 100:00:00
#SBATCH -J base_recalibration

module load bioinfo-tools
module load GATK/4.1.1.0
module load java/sun_jdk1.8.0_151
module load samtools


cd /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller

mkdir -p round2


for i in 101C_AGTTCC 102C_ATGTCA 1C_GCCAAT 2C_CAGATC 31C_ACTTGA 32C_GATCAG 61C_TAGCTT 62C_GGCTAC 91C_CTTGTA 92C_AGTCAA

do

  BAM=LsinapisSweM-Swe-sin-"$i"-sorted.deduped
  reference=LsinapisSweM.fasta

  cp /home/larshook/LarsH/FastZ/Snps/Read_mapping/MarkDuplicates/$BAM.bam $SNIC_TMP
  cp /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/$reference $SNIC_TMP
  cp /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/LsinapisSweM-Swe-sin_final_snps-filtered.vcf $SNIC_TMP


  cd $SNIC_TMP

  samtools index $BAM.bam
  samtools faidx $reference
  gatk --java-options "-Xmx6g" CreateSequenceDictionary -R $reference 
  gatk --java-options "-Xmx6g" IndexFeatureFile -F LsinapisSweM-Swe-sin_final_snps-filtered.vcf  



  gatk --java-options "-Xmx6g" BaseRecalibrator \
	-I $BAM.bam \
	-O $BAM.table \
	--known-sites LsinapisSweM-Swe-sin_final_snps-filtered.vcf \
	-R $reference

  gatk --java-options "-Xmx6g" ApplyBQSR \
	-I $BAM.bam \
	-O $BAM.bqsr.bam \
	--bqsr-recal-file $BAM.table \
	-R $reference

  cp $BAM.bqsr.bam /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2

done



for i in 1C_CCGTCC 21C_GTGAAA 22C_GTGGCC 2C_GTCCGC 41C_GTTTCG 42C_CGTACG 61C_GAGTGG 62C_ACTGAT 81C_ATTCCT 82C_ATCACG

do

  BAM=LsinapisSweM-Ire-juv-"$i"-sorted.deduped
  reference=LsinapisSweM.fasta

  cp /home/larshook/LarsH/FastZ/Snps/Read_mapping/MarkDuplicates/$BAM.bam $SNIC_TMP
  cp /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/$reference $SNIC_TMP
  cp /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/LsinapisSweM-Ire-juv_final_snps-filtered.vcf $SNIC_TMP


  cd $SNIC_TMP

  samtools index $BAM.bam
  samtools faidx $reference
  gatk --java-options "-Xmx6g" CreateSequenceDictionary -R $reference
  gatk --java-options "-Xmx6g" IndexFeatureFile -F LsinapisSweM-Ire-juv_final_snps-filtered.vcf

  gatk --java-options "-Xmx6g" BaseRecalibrator \
        -I $BAM.bam \
        -O $BAM.table \
        --known-sites LsinapisSweM-Ire-juv_final_snps-filtered.vcf \
        -R $reference

  gatk --java-options "-Xmx6g" ApplyBQSR \
        -I $BAM.bam \
        -O $BAM.bqsr.bam \
        --bqsr-recal-file $BAM.table \
        -R $reference

  cp $BAM.bqsr.bam /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2

done

