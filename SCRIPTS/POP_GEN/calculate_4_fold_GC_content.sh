#!/bin/bash -l

# calculate GC% at 4-fold degenerate sites from degenotate output, ignore N sites...

cd /home/larshook/LarsH/FastZ/PI

awk '{if ($5=="4") print}' degeneracy-all-sites_CDS.bed |\
	sed 's/-RA/ /' |\
		awk '{if ($7=="G" || $7=="C") {gc[$4]++;} else if ($7=="A" || $7=="T") {at[$4]++;}} END {for (i in gc) print i" "gc[i]/(gc[i]+at[i])}' |\
			sort -V -k 1,1 > LsinapisSweM_4_fold_gc.txt
