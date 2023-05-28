#!/bin/bash


for chr in {1..29}
do
  sbatch genotype_gvcf_second.sh $chr
  sleep 5
done
