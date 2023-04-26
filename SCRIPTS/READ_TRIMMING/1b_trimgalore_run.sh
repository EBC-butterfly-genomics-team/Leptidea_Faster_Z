#!bin/bash


for sample in 101 102 103 104 105 106 107 108
do 
  sbatch trimgalore_query.sh $sample
  sleep 5
done


cd /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS
mkdir -p 1b_TrimGalore
cd 1b_TrimGalore
mkdir -p REPORTS
