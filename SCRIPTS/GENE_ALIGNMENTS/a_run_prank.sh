#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 50:00:00
#SBATCH -J prank

module load bioinfo-tools
module load prank

main_path=/home/larshook/LarsH/FastZ/DN_DS

cd $main_path
mkdir -p ALIGNMENTS

for i in {00000001..00014378}

do
  prank \
	-d="$main_path"/ALIGNMENT_FASTA/gene_"$i".fasta \
	-o="$main_path"/ALIGNMENTS/gene_"$i" \
	-codon \
	-F \
	-iterate=1 \
	-tree="((LSSWET"$i", ((102_S23_"$i", 108_S29_"$i"), 106_S27_"$i")), 104_S25_"$i");"
done

# 102 - amurensis
# 104 - duponcheli
# 106 - lactea
# 108 - morsei
