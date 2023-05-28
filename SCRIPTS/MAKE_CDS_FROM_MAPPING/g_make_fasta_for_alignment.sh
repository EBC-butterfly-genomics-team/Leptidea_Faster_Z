#!/bin/bash -l


# concatenate genes from each species into one multifasta per gene to be used for alignment...


cd /home/larshook/LarsH/FastZ/DN_DS/
mkdir -p ALIGNMENT_FASTA

cd ORFS

# concatenate genes...

for i in {00000001..00014378}

do

  cat LSSWET"$i".fasta > ../ALIGNMENT_FASTA/gene_"$i".fasta
  cat 102_S23_"$i".fasta >> ../ALIGNMENT_FASTA/gene_"$i".fasta
  cat 104_S25_"$i".fasta >> ../ALIGNMENT_FASTA/gene_"$i".fasta
  cat 106_S27_"$i".fasta >> ../ALIGNMENT_FASTA/gene_"$i".fasta
  cat 108_S29_"$i".fasta >> ../ALIGNMENT_FASTA/gene_"$i".fasta


# check if all 5 species are present, if not remove fasta...

  gene_count=$(grep ">" ../ALIGNMENT_FASTA/gene_"$i".fasta | wc -l)
  
  if [ $gene_count != "5" ]
  then
    rm -f ../ALIGNMENT_FASTA/gene_"$i".fasta
  fi

done
