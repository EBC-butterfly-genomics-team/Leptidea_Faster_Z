#!/bin/bash


cd /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/SNPS
mkdir -p MPILEUP

cd /home/larshook/SCRIPTS/FastZ/Z_CHR_OUTGROUPS/Z_MF_SNP_DENSITY/MPILEUP


for species in P14458_103_S24 P14458_104_S25 P14458_107_S28 P14458_108_S29 P14502_103 P14502_104
do
  for chr in {1..29}
  do
    sbatch mpileup_query.sh $species $chr
    sleep 2
  done
  sleep 2
done
