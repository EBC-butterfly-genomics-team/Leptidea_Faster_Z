#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J recalculate_gene_coords

module load bioinfo-tools
module load emboss


# if a gene gets trimmed by orf prediction, recalculate the coordinates for downstream analyses...


cd /home/larshook/LarsH/FastZ/DN_DS

sp_id=LSSWET
seq_id=LSSWET
gff=/proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/gff/LsinapisSweM.all.maker.genes.gff
out_file=recalculated_gene_coordinates.txt


rm -f "$out_file"

for i in {00000001..00014378}
do

  # align trimmed orfs and untrimmed cds...

  needle \
        -asequence ORFS/"$sp_id""$i".fasta \
        -bsequence CDS/"$sp_id""$i".fasta \
        -gapopen 100 \
        -gapextend 10 \
        -outfile check_"$i".out \
        -aformat markx1



  gene_number="$seq_id""$i"

  # check for difference in start and end...
  
  # count gaps up until first start codon (A(TG))
  start_offset=$(grep -v "#" check_"$i".out | grep "$seq_id" | awk 'NR%2==1 {print $2}' | fold -w 1 | awk '{print $0};/A/{exit}' | wc -l | awk '{print $0-1}')
  
  # count gaps from end of sequence, skip if 3 gaps as these are removed stop codons (no need to recalculate coordinates)...
  end_offset=$(grep -v "#" check_"$i".out | grep "$seq_id" | awk 'NR%2==1 {print $2}' | fold -w 1 | tac | awk '{print $0};/[A-Z]/{exit}' | wc -l | awk '{if ($0==4) print 0; else print $0-1}')
  

  # if difference, recalculate coordinates...
  if [ $start_offset != "0" ] || [ $end_offset != "0" ]
  then
    
    # get chromosome...
    chr_id=$(grep -w "$gene_number" "$gff" | awk '{print $1}' | uniq)

    # make array with chr position and cds coordinate range...
    coord_array=( $(grep -w "$gene_number" "$gff" |\
	awk '{if ($3=="CDS") print $4".."$5}' |\
	while read line
	do 
	  for coord in $(eval echo {$line})
	  do 
	    echo $coord $line
	  done
	done |\
	sort -V -k 1,1 | awk '{print $1","$2}') )


    # count positions, and which positions to keep...
    cds_length=$(for position in "${coord_array[@]}"; do echo "$position"; done | wc -l)
    keep_start=$(echo $cds_length $start_offset | awk '{print $1-$2}')
    keep_end=$(echo $keep_start $end_offset | awk '{print $1-$2}')
    
    #keep_end=$(echo $cds_length $end_offset | awk '{print $1-$2}')


    # check orientation...
    gene_ori=$(grep -w "$gene_number" "$gff" | awk '{if ($3=="CDS") print $7}' | uniq )


    # trim positions...
    if [[ $gene_ori = "+" ]]
    then
      trimmed_coords=( $(for position in "${coord_array[@]}"; do echo "$position"; done |\
      tail -n "$keep_start" |\
      head -n "$keep_end") )
    
    elif [[ $gene_ori = "-" ]]
    then
      trimmed_coords=( $(for position in "${coord_array[@]}"; do echo "$position"; done |\
      head -n "$keep_start" |\
      tail -n "$keep_end") )
    fi

    
    # print new coords for each exon...
    for position in "${trimmed_coords[@]}"
    do
      echo $position
    done |\
      sed 's/,/ /' |\
      awk '{if (NR==1) print $1; else if ($2!=prev_exon) print prev_coord"\n"$1} {prev_exon=$2} {prev_coord=$1} END {print $1}' |\
      paste - - |\
      awk -v OFS="\t" -v chr_id=$chr_id -v ori=$gene_ori -v gene_number=$gene_number '{print chr_id, "recalculated", "CDS", $1, $2, ".", ori, "0", "Parent="gene_number"-RA"}' >> "$out_file"
  fi

  rm -f check_"$i".out

done
