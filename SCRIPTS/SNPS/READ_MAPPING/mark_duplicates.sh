#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 5:00:00
#SBATCH -J mark_duplicates

module load bioinfo-tools
module load GATK/4.1.1.0
module load java/sun_jdk1.8.0_151


BAM=LsinapisSweM-"$1"-sorted
FASTA=LsinapisSweM.fasta

cp /home/larshook/LarsH/FastZ/Snps/Read_mapping/$BAM.bam $SNIC_TMP
cp /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/$FASTA $SNIC_TMP

cd $SNIC_TMP

gatk --java-options "-Xmx12g" MarkDuplicates \
	-I $BAM.bam \
	-O $BAM.deduped.bam \
	-M $BAM.dup_metrics.txt \
	--ASSUME_SORT_ORDER coordinate \
	--CREATE_INDEX true \
	--REFERENCE_SEQUENCE $FASTA \
	--REMOVE_DUPLICATES true \
	--VERBOSITY INFO

cp $BAM.deduped.bam $BAM.dup_metrics.txt /home/larshook/LarsH/FastZ/Snps/Read_mapping/MarkDuplicates
