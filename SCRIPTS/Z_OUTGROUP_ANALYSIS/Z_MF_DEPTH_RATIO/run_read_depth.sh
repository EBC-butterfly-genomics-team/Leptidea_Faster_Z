#!/bin/bash


for i in LsinapisSweM-P14458_103_S24 LsinapisSweM-P14458_104_S25 LsinapisSweM-P14502_103 LsinapisSweM-P14502_104 LsinapisSweM-P14458_107_S28 LsinapisSweM-P14458_108_S29 LsinapisSweM-P14502_105 LsinapisSweM-P14502_106
do
  sbatch read_depth.sh $i
  sleep 2
done
