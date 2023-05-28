#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 15:00:00
#SBATCH -J filter_alignment


module load bioinfo-tools
module load hmmer
module load BioPerl


# filter alignments using HmmCleaner

msa_path=/home/larshook/LarsH/FastZ/DN_DS

cd "$msa_path"
mkdir -p FILTERED_ALIGNMENTS

for i in {00000001..00014378}

do

  /home/larshook/LarsH/SOFTWARE/Bio-MUST-Apps-HmmCleaner-0.180750/bin/HmmCleaner.pl "$msa_path"/ALIGNMENTS/gene_"$i".best.aa

  mv "$msa_path"/ALIGNMENTS/gene_"$i".best_hmm.* "$msa_path"/FILTERED_ALIGNMENTS

done
