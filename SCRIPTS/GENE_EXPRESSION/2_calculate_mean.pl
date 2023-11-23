#!/usr/bin/perl -w

use warnings; 
use strict;

# Lars Höök 2017
# Sum across samples and calculate mean FPKM for each gene


my ($outfile, $group_name) = @ARGV;

my %gene;

	open (OUT, ">>$outfile");
	print OUT "gene_id\tgene_name\tscaffold\tFPKM_$group_name\n";	#print headers

while(my $line = <STDIN>) {

	chomp $line;
	my @col = split (/\t/, $line); 					
	my $FPKM = pop @col;
	
	if($col[0] eq "gene_id") {next}

	my $gene = join '	' ,@col; 				#make hash key from all columns except FPKM
	$gene{$gene}{sample}++;						#sample count
	$gene{$gene}{sum} += $FPKM;					#sum FPKM's for same gene id

}


	
foreach(sort keys(%gene) ) {						#sort

	my $mean = $gene{$_}{sum}/$gene{$_}{sample};			#calculate mean

	open (OUT, ">>$outfile");   	
	print OUT "$_\t$mean\n";

}
