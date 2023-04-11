#!bin/bash


cd /home/larshook/LarsH/FastZ/
mkdir -p Snps
cd Snps
mkdir -p Read_trimming
cd Read_trimming
mkdir -p REPORTS


for sample in 101C 102C 1C 2C 31C 32C 61C 62C 91C 92C

do 
  sbatch trimgalore_query.sh $sample
  sleep 5
done
