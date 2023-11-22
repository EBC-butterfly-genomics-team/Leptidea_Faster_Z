#!/bin/bash -l
#SBATCH -J fastQC
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -A snic2022-5-34

# load modules
module load bioinfo-tools
module load FastQC/0.11.5
module load MultiQC/0.9

SRCDIR=/home/larshook/LarsH/RNAseq_Lsin_2020/1z_QC_concatenated_files

cd $SRCDIR


fastqc -o $SRCDIR/ *.gz

#fastqc -o $SRCDIR/ adjusted-condetri-P5052_202_S48_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_210_S51_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_218_S55_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_226_S59_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_202_S48_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_210_S51_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_218_S55_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_226_S59_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_203_S49_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_211_S52_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_219_S56_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_227_S38_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_203_S49_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_211_S52_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_219_S56_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_227_S38_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_204_S34_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_212_S36_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_220_S57_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P7553_325_S10_L002_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_204_S34_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_212_S36_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_220_S57_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P7553_325_S10_L002_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_205_S35_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_213_S53_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_221_S58_trim1.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_205_S35_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_213_S53_trim2.tagged_filter.fq.gz
#fastqc -o $SRCDIR/ adjusted-condetri-P5052_221_S58_trim2.tagged_filter.fq.gz


multiqc $SRCDIR/                           # creates single report


#note: runtime when using project folder (for 1 fasta.gz file): 6m
#                                        (for 12 fasta.gz.files, 3 to 4Gb each): 1h15m
#      runtime when using local scratch disk associated to node (for 1 fasta.gz file): 7m

