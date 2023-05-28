#!/bin/bash -l


# remove alignments if > 50% of initial alignment is trimmed...



main_path=/home/larshook/LarsH/FastZ/DN_DS

cd "$main_path"

mkdir -p TRUNCATED_ALIGNMENTS


for i in {00000001..00014378}
do

  # get sizes pre and post alignment...

  pre_size=$(fold -w 1 "$main_path"/ALIGNMENTS/gene_"$i".best.fas | wc -l)
  echo "$pre_size"

  post_size=$(fold -w 1 "$main_path"/TRUNCATED_ALIGNMENTS/gene_"$i".best.filtered.fas | wc -l)

  echo "$post_size"

  pp_frac=$(echo "scale=2; "$post_size" / "$pre_size"" | bc )

  # remove alignment if < 50%

  echo "$pp_frac"

  remove_align=$(echo "$pp_frac" | awk '{if ($1<0.50) print "yes"}')

  if [[ $remove_align = "yes" ]] 
  then
  rm -f "$main_path"/TRUNCATED_ALIGNMENTS/gene_"$i".best.filtered.fas
  fi

done
