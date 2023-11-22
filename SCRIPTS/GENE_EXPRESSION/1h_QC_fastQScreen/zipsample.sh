#!/bin/bash -l
#SBATCH -J zip
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -A snic2022-5-34


cp /home/larshook/LarsH/RNAseq_Lsin_2020/1x_QC_FINAL/adjusted-condetri-P7553_325_S10_L002_trim1.tagged_filter.fastq $SNIC_TMP
cp /home/larshook/LarsH/RNAseq_Lsin_2020/1x_QC_FINAL/adjusted-condetri-P7553_325_S10_L002_trim2.tagged_filter.fastq $SNIC_TMP

cd $SNIC_TMP

gzip -c $SNIC_TMP/adjusted-condetri-P7553_325_S10_L002_trim1.tagged_filter.fastq > adjusted-condetri-P7553_325_S10_L002_trim1.tagged_filter.fq.gz
gzip -c $SNIC_TMP/adjusted-condetri-P7553_325_S10_L002_trim2.tagged_filter.fastq > adjusted-condetri-P7553_325_S10_L002_trim2.tagged_filter.fq.gz

cp *.gz /home/larshook/LarsH/RNAseq_Lsin_2020/1z_QC_concatenated_files/
