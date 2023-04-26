#!/bin/bash -l

#SBATCH -J TrimGalore
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 35:00:00
#SBATCH -A snic2022-5-34


module load bioinfo-tools
module load TrimGalore/0.6.1
module load cutadapt/4.0
module load FastQC/0.11.9


cp /proj/uppstore2017185/b2014034/private/raw_data/Population_resequencing_data/151201_D00457_0128_AC80HPANXX/Sample_"$1"-"$2"/"$1"-"$2"*.fastq.gz $SNIC_TMP

cd $SNIC_TMP

trim_galore \
	--paired \
	--illumina \
	--phred33 \
	--stringency 1 \
	-e 0.1 \
	--nextseq 30 \
	--length 30 \
	--gzip \
	--fastqc \
	*R1_001.fastq.gz \
	*R2_001.fastq.gz

cp *fq.gz /home/larshook/LarsH/FastZ/Snps/Read_trimming/TRIMMED_READS
cp *txt *html /home/larshook/LarsH/FastZ/Snps/Read_trimming/TRIMMED_READS/REPORTS
