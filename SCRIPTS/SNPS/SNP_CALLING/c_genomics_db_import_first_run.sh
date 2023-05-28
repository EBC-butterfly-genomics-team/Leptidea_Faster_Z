#!/bin/bash

for chr in {1..29}
do
  sbatch genomics_db_import_first.sh $chr
  sleep 2
done
