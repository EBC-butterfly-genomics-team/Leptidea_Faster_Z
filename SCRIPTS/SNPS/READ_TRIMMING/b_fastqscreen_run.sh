#!/bin/bash -l



cd /home/larshook/LarsH/FastZ/Snps/Read_trimming/TRIMMED_READS

ls *.gz > /home/larshook/SCRIPTS/FastZ/SNPS/READ_TRIMMING/sample_list.txt

cd /home/larshook/SCRIPTS/FastZ/SNPS/READ_TRIMMING

while read line
do
  sbatch fastqscreen_query.sh $line
  sleep 5
done <sample_list.txt

cd /home/larshook/LarsH/FastZ/Snps/Read_trimming
mkdir -p FILTERED_READS
cd FILTERED_READS
mkdir -p REPORTS

rm -f /home/larshook/SCRIPTS/FastZ/SNPS/READ_TRIMMING/sample_list.txt
