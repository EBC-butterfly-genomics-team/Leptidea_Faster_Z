#!/bin/bash


for fold in 0 4 N
do
  sbatch pixy_fst.sh $fold
  sleep 2
done
