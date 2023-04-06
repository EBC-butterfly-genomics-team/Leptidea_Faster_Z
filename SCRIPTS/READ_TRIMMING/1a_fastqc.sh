#!/bin/bash -l

#SBATCH -A snic2022-5-34
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J fastQC

module load bioinfo-tools
module load FastQC/0.11.9
module load MultiQC/1.12

cd /home/larshook/LarsH/FastZ
mkdir -p LEPTIDEA_OUTGROUPS
cd LEPTIDEA_OUTGROUPS
mkdir -p 1a_fastqc_output
cd 1a_fastqc_output
outdir=$(pwd)

fastqc -o $outdir/ /proj/uppstore2017185/b2014034/private/raw_data/Population_resequencing_data/Leptidea_outgroups/P14458/P14458_101/02-FASTQ/191108_A00621_0143_BHTTVTDSXX/*.fastq.gz

fastqc -o $outdir/ /proj/uppstore2017185/b2014034/private/raw_data/Population_resequencing_data/Leptidea_outgroups/P14458/P14458_102/02-FASTQ/191108_A00621_0143_BHTTVTDSXX/*.fastq.gz

fastqc -o $outdir/ /proj/uppstore2017185/b2014034/private/raw_data/Population_resequencing_data/Leptidea_outgroups/P14458/P14458_103/02-FASTQ/191108_A00621_0143_BHTTVTDSXX/*.fastq.gz

fastqc -o $outdir/ /proj/uppstore2017185/b2014034/private/raw_data/Population_resequencing_data/Leptidea_outgroups/P14458/P14458_104/02-FASTQ/191108_A00621_0143_BHTTVTDSXX/*.fastq.gz

fastqc -o $outdir/ /proj/uppstore2017185/b2014034/private/raw_data/Population_resequencing_data/Leptidea_outgroups/P14458/P14458_105/02-FASTQ/191108_A00621_0143_BHTTVTDSXX/*.fastq.gz

fastqc -o $outdir/ /proj/uppstore2017185/b2014034/private/raw_data/Population_resequencing_data/Leptidea_outgroups/P14458/P14458_106/02-FASTQ/191108_A00621_0143_BHTTVTDSXX/*.fastq.gz

fastqc -o $outdir/ /proj/uppstore2017185/b2014034/private/raw_data/Population_resequencing_data/Leptidea_outgroups/P14458/P14458_107/02-FASTQ/191108_A00621_0143_BHTTVTDSXX/*.fastq.gz

fastqc -o $outdir/ /proj/uppstore2017185/b2014034/private/raw_data/Population_resequencing_data/Leptidea_outgroups/P14458/P14458_108/02-FASTQ/191108_A00621_0143_BHTTVTDSXX/*.fastq.gz


multiqc $outdir/
