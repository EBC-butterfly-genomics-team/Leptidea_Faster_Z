#!/bin/bash


printf "group1\tgroup2\tn1\tn2\tstatistic\tp\tp.adj\tp.adj.signif\n" dunn_test_dnds.txt > dunn_test_table.txt

sed 's/"//g;s/ //g' dunn_test_dnds.txt | awk -v OFS="\t" 'NR>1 {print $3, $4, $5, $6, $7, $8, $9, $10 }' >> dunn_test_table.txt


printf "chrN\tcount\n" > sign_dunn_test_count_for_plot.txt

grep "*" dunn_test_table.txt |\
        awk -F "\t" '{print $1"\n"$2}' |\
                awk -F "\t" '{c[$1]++} END {for (i in c) print i"\t"c[i]}' |\
                        sed 's/Chr/Chr /' >> sign_dunn_test_count_for_plot.txt
