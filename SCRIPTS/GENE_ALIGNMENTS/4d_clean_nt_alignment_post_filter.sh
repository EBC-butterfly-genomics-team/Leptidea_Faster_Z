#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 15:00:00
#SBATCH -J post_filter

module load java

msa_path=/home/larshook/LarsH/FastZ/DN_DS

for i in {00000001..00014378}
do
  java -jar /home/larshook/LarsH/SOFTWARE/macse_v2.06.jar -prog reportMaskAA2NT \
	-align "$msa_path"/ALIGNMENTS/gene_"$i".best.fas \
	-align_AA "$msa_path"/FILTERED_ALIGNMENTS/gene_"$i".best_hmm.fasta \
	-mask_AA - \
	-min_seq_to_keep_site 3 \
	-out_NT "$msa_path"/FILTERED_ALIGNMENTS/gene_"$i".best.filtered.fas

done
