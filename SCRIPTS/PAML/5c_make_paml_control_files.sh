#!/bin/bash

while read first second
do
  sed "s/xx/$first/;s/ndata = yy/ndata = $second/" block_template.ctl > genes_"$first"_block.ctl
done < /home/larshook/LarsH/FastZ/DN_DS/FILTERED_ALIGNMENTS/number_of_genes_per_block.txt
