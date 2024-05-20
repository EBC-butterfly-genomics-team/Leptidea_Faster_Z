#!/bin/bash -l

#SBATCH -A naiss2024-5-14
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 10:00:00
#SBATCH -J paml


module load bioinfo-tools
module load paml/4.9j


cp /home/larshook/LarsH/FastZ/DN_DS/FILTERED_ALIGNMENTS/genes_"$1"_block.paml $SNIC_TMP
cp /home/larshook/SCRIPTS/FastZ/PAML/genes_"$1"_block.ctl $SNIC_TMP 
cp /home/larshook/SCRIPTS/FastZ/PAML/branch_lep_tree.txt $SNIC_TMP

cd $SNIC_TMP

codeml genes_"$1"_block.ctl

cp *.out /home/larshook/LarsH/FastZ/DN_DS/MULTI_BRANCH/
