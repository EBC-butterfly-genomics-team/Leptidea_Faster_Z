#!/bin/bash -l
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 15:00:00
#SBATCH -J genotype_gvcf

module load bioinfo-tools
module load GATK/4.1.1.0
module load java/sun_jdk1.8.0_151
module load samtools

reference=LsinapisSweM
vcf_path=/home/larshook/LarsH/FastZ/Snps/HaplotypeCaller

cp /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/$reference.fasta $SNIC_TMP

cd $SNIC_TMP

samtools faidx $reference.fasta
gatk --java-options "-Xmx6g" CreateSequenceDictionary -R $reference.fasta

for species in Swe-sin Ire-juv
do
  for chr in {1..29}
  do

    cp -r $vcf_path/"$reference"-"$species"_HiC_scaffold_"$chr" $SNIC_TMP

    gatk --java-options "-Xmx6g" GenotypeGVCFs \
	-R $reference.fasta \
	-V gendb://$reference-"$species"_HiC_scaffold_"$chr" \
	-O $reference-"$species"_HiC_scaffold_"$chr".vcf

    cp $reference-"$species"_HiC_scaffold_"$chr"* $vcf_path

  done
done
