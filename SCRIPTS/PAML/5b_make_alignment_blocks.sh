#!/bin/bash


# make blocks of gene alignments for input to paml

genes_per_block=2500


cd /home/larshook/LarsH/FastZ/DN_DS/FILTERED_ALIGNMENTS


ls *filtered.phy > gene_list.txt

split -l "$genes_per_block" --numeric-suffixes gene_list.txt genes_


wc -l genes_* | awk '{print $2, $1}' | sed 's/genes_//' | grep -v "total" > number_of_genes_per_block.txt


for file in genes_*
do
  while read line
  do
    cat "$line" >> "$file"_block.paml
  done <"$file"
done
