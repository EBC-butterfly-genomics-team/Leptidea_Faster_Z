#!/usr/bin/perl -w

use warnings; 
use strict;

# Lars Höök 2017
# Sort by gene id, sum FPKM's for each gene id (if more than one), and remove unused columns


my ($infile, $outfile) = @ARGV ;

my %gene;

open(IN, $infile) or die "Can't open $infile\n";

while(my $line = <IN>) {
	chomp $line;
	my @col = split (/\t/, $line); 					
	$col[1] =~ s/ (?<!\d) \. (?!\d) /-/xg;				#replace . with - for concatenation
	pop @col;							#remove column: "TPM"
	my $FPKM = pop @col;
	splice @col, 3, 4; 						#remove columns: "Strand", "Start", "End", "Coverage"

	if($col[0] eq "Gene_ID") {

		open (OUT, ">>$outfile");
		print OUT "gene_id\tgene_name\tscaffold\tFPKM\n";	#print headers
	
		next

	}

	my $gene = join '	' ,@col; 				#make hash key from all columns except FPKM
	$gene{$gene}{vol} += $FPKM;					#sum FPKM's for same gene id

}

	
foreach(sort keys(%gene) ) {						#sort and print

	open (OUT, ">>$outfile");   	
	print OUT "$_\t$gene{$_}{vol}\n";

}
