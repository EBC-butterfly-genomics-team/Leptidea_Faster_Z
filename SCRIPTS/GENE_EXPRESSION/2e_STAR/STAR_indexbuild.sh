#!/bin/bash -l

#SBATCH -J Indexbuild
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 1:00:00
#SBATCH -A naiss2023-5-52

# load modules
module load bioinfo-tools
module load star/2.7.9a


STAR \
    --runMode genomeGenerate \
    --runThreadN 8 \
    --genomeDir /home/larshook/LarsH/FastZ/RNAseq_Lsin/STAR_index \
    --genomeFastaFiles /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/LsinapisSweM.fasta \
    --genomeSAindexNbases 13
