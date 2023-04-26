#!/bin/bash -l

#SBATCH -A snic2022-5-34
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J BBmap_repair

module load bioinfo-tools
module load bbmap/38.61b
module load gzip

cd /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS
mkdir -p READS_FINAL

for i in 101_S22 102_S23 103_S24 104_S25 105_S26 106_S27 107_S28 108_S29
do 

  cp /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/FILTERED_READS/P14458_"$i"_L003_R*_001_val_*.tagged_filter.fastq.gz $SNIC_TMP

  cd $SNIC_TMP

  repair.sh -Xmx6G \
        in1=P14458_"$i"_L003_R1_001_val_1.tagged_filter.fastq.gz \
	in2=P14458_"$i"_L003_R2_001_val_2.tagged_filter.fastq.gz \
	out1=P14458_"$i"_L003_R1_001_val_1_final.fq \
	out2=P14458_"$i"_L003_R2_001_val_2_final.fq \
	outs=P14458_"$i"_L003_singeltons_final.fq \
	repair


  gzip *final.fq

  cp *final.fq.gz /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READS_FINAL

done
