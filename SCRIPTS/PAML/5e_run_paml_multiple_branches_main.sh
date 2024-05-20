#!/bin/bash

#for i in {00..11}
#for i in 01 02 05
for i in 01
do
  sbatch paml_mb_query.sh "$i"
  sleep 2
done

#cd /home/larshook/LarsH/FastZ/DN_DS/
#mkdir -p MULTI_BRANCH
