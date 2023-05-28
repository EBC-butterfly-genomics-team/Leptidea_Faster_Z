#!/bin/bash


module load bioinfo-tools
module load BEDTools


cd /home/larshook/LarsH/FastZ


# make 100kb windows...
bedtools makewindows -g LsinapisSweM.genome -w 100000 > LsinapisSweM_100kb.bed
