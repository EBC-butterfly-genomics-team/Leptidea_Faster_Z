#!/bin/bash

for i in LsinapisSweF-P14458_103_S24 LsinapisSweF-P14458_104_S25 LsinapisSweF-P14502_103 LsinapisSweF-P14502_104 LsinapisSweF-P14458_107_S28 LsinapisSweF-P14458_108_S29
do
  sbatch read_depth.sh $i
  sleep 5
done
