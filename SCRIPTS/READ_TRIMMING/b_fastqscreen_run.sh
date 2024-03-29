#!/bin/bash -l


cd /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/TRIMMED_READS

ls *.gz > /home/larshook/SCRIPTS/FastZ/READ_TRIMMING/sample_list.txt

cd /home/larshook/SCRIPTS/FastZ/READ_TRIMMING

while read line
do
  sbatch fastqscreen_query.sh $line
  sleep 5
done <sample_list.txt

cd /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS
mkdir -p FILTERED_READS
cd FILTERED_READS
mkdir -p REPORTS

rm -f /home/larshook/SCRIPTS/FastZ/READ_TRIMMING/sample_list.txt
