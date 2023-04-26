#!bin/bash


for sample in 101 102 103 104 105 106 107 108
do 
  sbatch trimgalore_query.sh $sample
  sleep 5
done


cd /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS
mkdir -p TRIMMED_READS
cd TRIMMED_READS
mkdir -p REPORTS
