#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J make_alignment_fasta


# concatenate genes from each species into one multifasta to be used for alignment...


module load bioinfo-tools
module load bbmap


cd /home/larshook/LarsH/FastZ/DN_DS
mkdir -p ALIGNMENT_FASTA
cd ORFS

# remove empty fasta

find . -size 0 -delete

# remove genes with more than 50% N's...

for file in *.fasta
do
  gene_n=$(stats.sh in="$file" | awk 'NR==2 {if ($5>0.5) print "skip_gene"}')

  if [[ $gene_n = "skip_gene" ]]
  then
    rm -f "$file"
  fi
done


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
