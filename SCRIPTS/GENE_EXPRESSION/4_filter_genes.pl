#!/usr/bin/perl -w

use strict;
use warnings;

#Lars Höök 2017
#Filter genes based on expression

#Genes are filtered out if they are below a defined treshold in any sex.
#The gene has to be expressed above defined threshold in both sexes to pass the filter.
#This removes genes with low and sex-specific expression.

############################# Apply filter here #######################################

my $filter = 0;
#my $filter = 1;

#######################################################################################


my ($infile, $prefix, $outfile) = @ARGV ;
my $logfile = "prepare_stringtie_files.log";
my $totalgenes;
my $filteredgenes;


open(IN, $infile) or die "Can't open $infile\n" ;

while (my $line = <IN>) {				#read through genes

	my @col = split(/\s+/, $line);
	
	if ($col[0] eq "gene_id") {			#print headers
		
	open (OUT, ">>$outfile") ;
	print OUT "$line" ;

	}

	else {
	
		if ($col[3] > $filter) {		#check first group
			
			if ($col[4] > $filter) {	#check second group

			open (OUT, ">>$outfile");
			print OUT "$line";
			$totalgenes++;			#log total genes

			}
			
			else {

			$filteredgenes++;		#log filtered genes
			next
			
			}			

		}		

		else {

		$filteredgenes++;			#log filtered genes
		next

		}

	}

}

close (IN);
close (OUT);

open (LOG, ">>$logfile");				#add to logfile
print LOG 
	"Genes $prefix filtered out: $filteredgenes\n",
	"Total genes $prefix for analysis: $totalgenes\n";

close (LOG);

