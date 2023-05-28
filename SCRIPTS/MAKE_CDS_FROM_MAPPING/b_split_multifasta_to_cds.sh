#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J split_fasta


# split the L. sinapis cds multifasta (ref, var, and snp masked) to separate cds before orf prediction...

# snp masked...

cp /home/larshook/LarsH/FastZ/LsinapisSweM_snp_masked_concatenated_CDS.fasta $SNIC_TMP

cd $SNIC_TMP

awk '{print $1}' LsinapisSweM_snp_masked_concatenated_CDS.fasta |\
   awk '/^>/ {OUT=substr($0,2) ".fasta"}; OUT {print >OUT}'

sleep 10m

# get chr genes only...

for i in {00000001..00014378}
do
  cp LSSWET"$i".fasta /home/larshook/LarsH/FastZ/DN_DS/CDS
done


# reference...

cp /home/larshook/LarsH/FastZ/LsinapisSweM_concatenated_CDS.fasta $SNIC_TMP

cd $SNIC_TMP

awk '{print $1}' LsinapisSweM_concatenated_CDS.fasta |\
   awk '/^>/ {OUT=substr($0,2) ".fasta"}; OUT {print >OUT}'

sleep 10m

# get chr genes only...

for i in {00000001..00014378}
do
  cp LSSWET"$i".fasta /home/larshook/LarsH/FastZ/PN_PS/CDS/LSREFT"$i".fasta
done


# variant...

cp /home/larshook/LarsH/FastZ/LsinapisSweM_variant_concatenated_CDS.fasta $SNIC_TMP

cd $SNIC_TMP

awk '{print $1}' LsinapisSweM_variant_concatenated_CDS.fasta |\
   awk '/^>/ {OUT=substr($0,2) ".fasta"}; OUT {print >OUT}'

sleep 10m

# get chr genes only...

for i in {00000001..00014378}
do
  cp LSSWET"$i".fasta /home/larshook/LarsH/FastZ/PN_PS/CDS/LSVART"$i".fasta
done
