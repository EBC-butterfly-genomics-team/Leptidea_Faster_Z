#!/bin/bash -l
#SBATCH -J catting
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -A snic2022-5-34

# load modules
module load bioinfo-tools

#input data files

PATH_MAIN=/home/larshook/LarsH/RNAseq_Lsin_2020       #project folder           
PATH_LOCAL=1x_QC_FINAL			             	             #folder containing input data

# output folder
OUT_FOLDER=1z_QC_concatenated_files

SRCDIR=$(pwd)                                                #remember current path

cd $PATH_MAIN/$PATH_LOCAL

cd $SRCDIR

# concatenate files associated to the same sample (applies to lane 5 and 6)

cp $PATH_MAIN/$PATH_LOCAL/adjusted-condetri-P5052_226_S59_L006_trim1.tagged_filter.fastq $SNIC_TMP
cp $PATH_MAIN/$PATH_LOCAL/adjusted-condetri-P5052_226_S59_L006_trim2.tagged_filter.fastq $SNIC_TMP
cp $PATH_MAIN/$PATH_LOCAL/adjusted-condetri-P5052_226_S59_L005_trim1.tagged_filter.fastq $SNIC_TMP
cp $PATH_MAIN/$PATH_LOCAL/adjusted-condetri-P5052_226_S59_L005_trim2.tagged_filter.fastq $SNIC_TMP

cat $SNIC_TMP/adjusted-condetri-P5052_226_S59_L005_trim1.tagged_filter.fastq $SNIC_TMP/adjusted-condetri-P5052_226_S59_L006_trim1.tagged_filter.fastq > $SNIC_TMP/adjusted-condetri-P5052_226_S59_trim1.tagged_filter.fastq

cat $SNIC_TMP/adjusted-condetri-P5052_226_S59_L005_trim2.tagged_filter.fastq $SNIC_TMP/adjusted-condetri-P5052_226_S59_L006_trim2.tagged_filter.fastq > $SNIC_TMP/adjusted-condetri-P5052_226_S59_trim2.tagged_filter.fastq

gzip -c $SNIC_TMP/adjusted-condetri-P5052_226_S59_trim1.tagged_filter.fastq > $PATH_MAIN/$OUT_FOLDER/adjusted-condetri-P5052_226_S59_trim1.tagged_filter.fq.gz
gzip -c $SNIC_TMP/adjusted-condetri-P5052_226_S59_trim2.tagged_filter.fastq > $PATH_MAIN/$OUT_FOLDER/adjusted-condetri-P5052_226_S59_trim2.tagged_filter.fq.gz





