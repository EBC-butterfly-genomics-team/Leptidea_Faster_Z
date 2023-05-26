#!/bin/bash

for sample in P14458_101_S22 P14458_102_S23 P14458_103_S24 P14458_104_S25 P14458_105_S26 P14458_106_S27 P14458_107_S28 P14458_108_S29
do
  sbatch bwa_query.sh $sample
  sleep 5
done

cd /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS
mkdir -p READ_MAPPING
