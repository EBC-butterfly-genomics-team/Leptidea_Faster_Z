#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J concatenate_cds

# Concatenate all cds from each gene into a single fasta sequence

module load bioinfo-tools
module load biopython
module load seqtk

cd /home/larshook/LarsH/FastZ
mkdir -p CDS_CONSENSUS

script_path=$(pwd)
CDS_path=/home/larshook/LarsH/FastZ/CDS_CONSENSUS
gene_id=LSSWET
reference=LsinapisSweM

cd $CDS_path
mkdir -p concatenated_cds

cd $script_path

for species_id in 102_S23 104_S25 106_S27 108_S29
#for species_id in 103_S24 107_S28

do

  # select only chr genes...
  for i in {00000001..00014378}
  #for i in {00000001..00016294}

  do

    # select one gene at a time and put cds in temporary fasta, using script "pick_scaffolds_from_multifasta.py"

    grep ""$gene_id""$i"" "$CDS_path"/"$reference"_CDS_regions_and_gene_info.txt |\
#   grep ""$gene_id""$i"" "$CDS_path"/recalculated_"$reference"_CDS_regions_and_gene_info.txt > "$CDS_path"/concatenated_cds/list"$i".txt
	awk -v species_id="$species_id" '{print species_id"_"$1}' > "$CDS_path"/concatenated_cds/list"$i".txt

    python pick_scaffolds_from_multifasta.py "$CDS_path"/"$species_id"_concensus_CDS.fa "$CDS_path"/concatenated_cds/list"$i".txt "$CDS_path"/concatenated_cds/temp_"$i".fasta


    # print fasta header based on gene name

    printf ">""$species_id""_""$i""\n" > "$CDS_path"/concatenated_cds/"$species_id"_"$i".fasta


    #check if gene is +/- strand, reverse complement if -

    str_ori=$(grep ""$gene_id""$i"" "$CDS_path"/"$reference"_CDS_regions_and_gene_info.txt |\
#   str_ori=$(grep ""$gene_id""$i"" "$CDS_path"/recalculated_"$reference"_CDS_regions_and_gene_info.txt |\
	awk '{print $2}' |\
		uniq)
	      	

    if [[ $str_ori = "+" ]]
    then
      grep -vE "^>" "$CDS_path"/concatenated_cds/temp_"$i".fasta |\
		awk '{if($0 ~ /^>/) {print "\n"$0} else {printf $0"\n"}}' >> "$CDS_path"/concatenated_cds/"$species_id"_"$i".fasta
    elif [[ $str_ori = "-" ]]
    then
      grep -vE "^>" "$CDS_path"/concatenated_cds/temp_"$i".fasta |\
		awk '{if($0 ~ /^>/) {print "\n"$0} else {printf $0"\n"}}' >> "$CDS_path"/concatenated_cds/"$species_id"_"$i".fasta
      seqtk seq -r "$CDS_path"/concatenated_cds/"$species_id"_"$i".fasta > "$CDS_path"/concatenated_cds/"$species_id"_"$i"_rev.fasta
      mv "$CDS_path"/concatenated_cds/"$species_id"_"$i"_rev.fasta "$CDS_path"/concatenated_cds/"$species_id"_"$i".fasta
    fi

    # remove temporary fasta and list

    rm -f "$CDS_path"/concatenated_cds/list"$i".txt
    rm -f "$CDS_path"/concatenated_cds/temp_"$i".fasta

  done
done
