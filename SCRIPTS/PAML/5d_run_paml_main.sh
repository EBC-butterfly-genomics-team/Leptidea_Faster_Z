#!/bin/bash

for i in {00..04}
do
  sbatch paml_query.sh "$i"
  sleep 2
done
