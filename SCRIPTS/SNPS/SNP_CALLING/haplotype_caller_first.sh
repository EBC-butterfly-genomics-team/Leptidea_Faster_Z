#!/bin/bash -l
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J haplotype_caller

module load bioinfo-tools
module load GATK/4.1.1.0
module load java/sun_jdk1.8.0_151
module load samtools


BAM=LsinapisSweM-"$1"-sorted.deduped
FASTA=LsinapisSweM.fasta

cp /home/larshook/LarsH/FastZ/Snps/Read_mapping/MarkDuplicates/$BAM.bam $SNIC_TMP
cp /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/$FASTA $SNIC_TMP

cd $SNIC_TMP
  
samtools index $BAM.bam
samtools faidx $FASTA

gatk --java-options "-Xmx6g" CreateSequenceDictionary -R $FASTA

gatk --java-options "-Xmx6g" HaplotypeCaller \
    -I $BAM.bam \
    -R $FASTA \
    -ERC GVCF \
    -L HiC_scaffold_"$2" \
    -O "$BAM"_HiC_scaffold_"$2".g.vcf.gz

cp "$BAM"_HiC_scaffold_"$2".g.vcf.gz /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller
