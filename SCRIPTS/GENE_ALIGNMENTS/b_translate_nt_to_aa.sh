#!/bin/bash -l

module load bioinfo-tools
module load java

cd /home/larshook/LarsH/FastZ/DN_DS/ALIGNMENTS

for i in {00000001..00014378}
do
  java -jar /home/larshook/LarsH/SOFTWARE/macse_v2.06.jar -prog translateNT2AA \
	-seq gene_"$i".best.fas \
	-out_AA gene_"$i".best.aa
done
