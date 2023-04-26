#!/bin/bash -l
#SBATCH -J fastqscreen
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 40:00:00
#SBATCH -A naiss2023-5-52


module load bioinfo-tools
module load fastq_screen/0.11.1
module load bwa/0.7.17


fastq_screen \
        --threads 5 \
        --aligner bwa \
	--conf fastq_screen.conf \
        --subset 0 \
        --nohits \
        --outdir /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/FILTERED_READS \
	$1

cd /home/larshook/LarsH/FastZ/Snps/Read_trimming/FILTERED_READS
mv *txt *png *html REPORTS/
rm -f *trimmed.tagged.fastq.gz
