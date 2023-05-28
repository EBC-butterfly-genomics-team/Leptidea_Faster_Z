#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J concatenate_cds

# Concatenate all cds from each gene into a single fasta sequence/file

module load bioinfo-tools
module load biopython
module load seqtk


script_path=$(pwd)
CDS_path=/home/larshook/LarsH/FastZ/DN_DS
gene_id=LSSWET
reference=LsinapisSweM
out_folder=CDS


cd $CDS_path
mkdir -p $out_folder


cd $script_path

for species_id in 102_S23 104_S25 106_S27 108_S29
do

  # select only chr genes...
  for i in {00000001..00014378}
  do

    # select one gene at a time and put cds in temporary fasta, using script "pick_scaffolds_from_multifasta.py"

    grep ""$gene_id""$i"" "$CDS_path"/"$reference"_CDS_regions_and_gene_info.txt |\
	awk -v species_id="$species_id" '{print species_id"_"$1}' > "$CDS_path"/$out_folder/list"$i".txt

    python pick_scaffolds_from_multifasta.py "$CDS_path"/"$species_id"_CDS.fa "$CDS_path"/$out_folder/list"$i".txt "$CDS_path"/$out_folder/temp_"$i".fasta


    # print fasta header based on gene name

    printf ">""$species_id""_""$i""\n" > "$CDS_path"/$out_folder/"$species_id"_"$i".fasta


    #check if gene is +/- strand, reverse complement if -

    str_ori=$(grep ""$gene_id""$i"" "$CDS_path"/"$reference"_CDS_regions_and_gene_info.txt |\
	awk '{print $2}' |\
		uniq)
	      	

    if [[ $str_ori = "+" ]]
    then
      grep -vE "^>" "$CDS_path"/$out_folder/temp_"$i".fasta |\
		awk '{if($0 ~ /^>/) {print "\n"$0} else {printf $0"\n"}}' >> "$CDS_path"/$out_folder/"$species_id"_"$i".fasta
    elif [[ $str_ori = "-" ]]
    then
      grep -vE "^>" "$CDS_path"/$out_folder/temp_"$i".fasta |\
		awk '{if($0 ~ /^>/) {print "\n"$0} else {printf $0"\n"}}' >> "$CDS_path"/$out_folder/"$species_id"_"$i".fasta
      seqtk seq -r "$CDS_path"/$out_folder/"$species_id"_"$i".fasta > "$CDS_path"/$out_folder/"$species_id"_"$i"_rev.fasta
      mv "$CDS_path"/$out_folder/"$species_id"_"$i"_rev.fasta "$CDS_path"/$out_folder/"$species_id"_"$i".fasta
    fi

    # remove temporary fasta and list

    rm -f "$CDS_path"/$out_folder/list"$i".txt
    rm -f "$CDS_path"/$out_folder/temp_"$i".fasta

  done
done



# remove empty fasta files

cd "$CDS_path"/$out_folder
find - size 18c -delete
find - size 19c	-delete
