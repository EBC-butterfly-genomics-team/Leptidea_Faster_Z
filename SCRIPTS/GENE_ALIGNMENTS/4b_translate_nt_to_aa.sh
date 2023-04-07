#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J translate_to_aa

module load bioinfo-tools
module load java

cd /home/larshook/LarsH/FastZ/DN_DS/ALIGNMENTS

for i in {00000001..00014378}
do
  java -jar /home/larshook/LarsH/SOFTWARE/macse_v2.06.jar -prog translateNT2AA \
	-seq gene_"$i".best.fas \
	-out_AA gene_"$i".best.aa
done
