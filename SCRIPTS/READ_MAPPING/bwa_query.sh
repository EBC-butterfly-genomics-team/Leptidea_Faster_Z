#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 150:00:00
#SBATCH -J bwa

module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.14

ASSEMBLY=LsinapisSweM
#ASSEMBLY=LsinapisSweF

lane=3

cp /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/"$ASSEMBLY".fasta $SNIC_TMP
cp /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READS_FINAL/"$1"_L003_R*_001_val_*_final.fq.gz $SNIC_TMP
#cp /home/larshook/LarsH/W_CHROMOSOME/READS_FINAL/P14502/"$1".barcoded_R*_final.fq.gz $SNIC_TMP

cd $SNIC_TMP

bwa index $ASSEMBLY.fasta

bwa mem \
	-M \
	-t 4 \
	-R "@RG\tID:"$1"\tSM:"$1"\tLB:$1\tPL:ILLUMINA\tPU:HTTVTDSXX_"$lane"" \
	$ASSEMBLY.fasta \
	"$1"_L003_R1_001_val_1_final.fq.gz \
	"$1"_L003_R2_001_val_2_final.fq.gz | samtools sort -@4 > $ASSEMBLY-$1-sorted.bam

#        "$1".barcoded_R1_final.fq.gz \
#        "$1".barcoded_R2_final.fq.gz | samtools sort -@4 > $ASSEMBLY-$1-sorted.bam


samtools index $ASSEMBLY-$1-sorted.bam $ASSEMBLY-$1-sorted.bai
samtools stats $ASSEMBLY-$1-sorted.bam > $ASSEMBLY-$1-sorted-stats.txt
samtools flagstat -O tsv $ASSEMBLY-$1-sorted.bam > $ASSEMBLY-$1-sorted-flagstats.txt

samtools view -h -@4 $ASSEMBLY-$R1-sorted.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b -q 30 -@4 > $ASSEMBLY-$1-sorted-unique.bam
samtools index $ASSEMBLY-$1-sorted-unique.bam $ASSEMBLY-$1-sorted-unique.bai

samtools stats $ASSEMBLY-$1-sorted-unique.bam > $ASSEMBLY-$1-sorted-unique-stats.txt
samtools flagstat -O tsv $ASSEMBLY-$1-sorted-unique.bam > $ASSEMBLY-$1-sorted-unique-flagstat.txt

cp *.bam *.bai *.txt /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_MAPPING
