#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 30:00:00
#SBATCH -J bwa

module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.14

ASSEMBLY=LsinapisSweM

cp /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/"$ASSEMBLY".fasta $SNIC_TMP
cp /home/larshook/LarsH/FastZ/Snps/Read_trimming/READS_FINAL/"$1"*R*final.fq.gz $SNIC_TMP


cd $SNIC_TMP

bwa index $ASSEMBLY.fasta

bwa mem \
	-M \
	-t 1 \
	-R "@RG\tID:"$1"\tSM:"$1"\tLB:$1\tPL:ILLUMINA\tPU:C80HPANXX_5" \
	$ASSEMBLY.fasta \
	"$1"*R1_final.fq.gz \
	"$1"*R2_final.fq.gz | samtools sort -@1 > $ASSEMBLY-$1-sorted.bam

samtools index $ASSEMBLY-$1-sorted.bam $ASSEMBLY-$1-sorted.bai
samtools stats $ASSEMBLY-$1-sorted.bam > $ASSEMBLY-$1-sorted-stats.txt
samtools flagstat -O tsv $ASSEMBLY-$1-sorted.bam > $ASSEMBLY-$1-sorted-flagstats.txt

cp *.bam *.bai *.txt /home/larshook/LarsH/FastZ/Snps/Read_mapping
