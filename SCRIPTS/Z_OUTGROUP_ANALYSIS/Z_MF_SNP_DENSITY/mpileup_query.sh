#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J mpileup

module load bioinfo-tools
module load bcftools
module load samtools

BAM=LsinapisSweM-"$1"-sorted.bam-unique.deduped
reference=LsinapisSweM.fasta

cp /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_MAPPING/MARK_DUPLICATES/$BAM*bam $SNIC_TMP
cp /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/$reference $SNIC_TMP


cd $SNIC_TMP

samtools index $BAM.bam
samtools faidx $reference


bcftools mpileup \
	-Ou \
	--count-orphans \
	--skip-indels \
	-r HiC_scaffold_"$2" \
	--fasta-ref $reference \
	$BAM.bam |\
		bcftools call \
			-mv \
			-Ob \
			-o $1-HiC_scaffold_"$2".bcf

cp $1-HiC_scaffold_"$2".bcf /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/SNPS/MPILEUP
