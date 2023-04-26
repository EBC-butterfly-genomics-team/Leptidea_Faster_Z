#!bin/bash


for sample in 101C 102C 1C 2C 31C 32C 61C 62C 91C 92C

do 
  sbatch trimgalore_query.sh Swe-sin $sample
  sleep 5
done


for sample in 1C 21C 22C 2C 41C 42C 61C 62C 81C 82C

do
  sbatch trimgalore_query.sh Ire-juv $sample
  sleep 5
done


cd /home/larshook/LarsH/FastZ/
mkdir -p Snps
cd Snps
mkdir -p Read_trimming
cd Read_trimming
mkdir -p REPORTS
