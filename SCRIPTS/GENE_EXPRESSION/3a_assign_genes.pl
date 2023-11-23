#!/usr/bin/perl -w

use warnings;
use strict;

# Lars Höök 2017
# Assign genes to A or Z

my ($scaffold_list, $gene_list) = @ARGV ;

my $prefix = $gene_list;
	$prefix =~ s/\-concatenated.txt//;

my $outfile = "$prefix-assigned_A_or_Z.txt";
my $logfile = "prepare_stringtie_files.log";

my $totalgenes;
my $assignedgenes;
my $sum_A;
my $sum_Z;

open (GENE, $gene_list) or die "Can't open $gene_list\n";
	
while (my $gene = <GENE>) {					#read through genes 

	chomp ($gene);
	my @col_2 = split(/\t/, $gene); 			
	my $gene_id = $col_2[0]; 
	my $match = $col_2[2];
		
	if ($gene_id eq "gene_id") {				#print headers
		
		open (OUT, ">>$outfile");
		print OUT "$gene\tchromosome\n";

	}
	
	else {$totalgenes++}					#log total genes

	open (IN, $scaffold_list) or die "Can't open $scaffold_list\n";

	while (my $line = <IN>) {				#read scaffolds

		my @col_1 = split(/\s+/, $line);		
		my $scaffold = shift(@col_1);
		my $chromosome = pop @col_1;
	
	
		if ($scaffold eq $match) {			#add gene to outfile if scaffold name matches
					
			$assignedgenes++;			#log assigned genes
			if ($chromosome eq "A") {$sum_A++}	#log A genes
			elsif ($chromosome eq "Z") {$sum_Z++}	#log Z genes

			open(OUT, ">>$outfile");
			print OUT "$gene\t$chromosome\n";	
			
		}			
		
		else { next }	

	}

}

close (IN);
close (GENE);
close (OUT);

my $nonassigned = $totalgenes-$assignedgenes;			#make log file
open (LOG, ">>$logfile");
print LOG "\n",
	"##############################\n",
	"\n",
	"Stage: $prefix\n", "\n",
	"Total genes: $totalgenes\n",
	"Assigned genes: $assignedgenes\n", 
	"Assigned A: $sum_A\n",
	"Assigned Z: $sum_Z\n",
	"Non-assigned genes: $nonassigned\n";

close (LOG);

