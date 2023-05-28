#!/bin/bash


# remove gene if not all 5 species are present after filtering


msa_path=/home/larshook/LarsH/FastZ/DN_DS/FILTERED_ALIGNMENTS


cd "$msa_path"

genes_before=$(ls gene_*.best.filtered.fas | wc -l)

for i in {00000001..00014378}
do

  # check if file exists...

  file=gene_"$i".best.filtered.fas

  if [ -f "$file" ];
  then

    # skip gene if not all species present...

    remove=$(grep ">" "$msa_path"/gene_"$i".best.filtered.fas | wc -l | awk '{if ($1!=5) print "yes"}')

    if [[ $remove = "yes" ]]
    then
      rm -f "$msa_path"/gene_"$i".best.filtered.fas
    fi
  fi

done

sleep 2

genes_after=$(ls gene_*.best.filtered.fas | wc -l)

printf "\ngenes before: $genes_before\n"
printf "genes after: $genes_after\n\n"
