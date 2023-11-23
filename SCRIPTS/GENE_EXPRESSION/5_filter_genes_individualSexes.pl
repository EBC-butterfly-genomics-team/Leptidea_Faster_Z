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


my ($infile, $prefix, $refsex, $outfile) = @ARGV ;
my $logfile = "prepare_stringtie_files.log";
my $totalgenes;
my $filteredgenes;


open(IN, $infile) or die "Can't open $infile\n" ;

if ($refsex eq "female") {				#read reference sex

        my $ref_col = 3

        }
        
        elsif ($refsex eq "male") {				

        my $ref_col = 4

        }




while (my $line = <IN>) {				#read through genes

	my @col = split(/\s+/, $line);
	
	if ($col[0] eq "gene_id") {			#print headers
		
	open (OUT, ">>$outfile") ;
	print OUT "$line" ;

	}

	else {
	
		if ($refsex eq "female") {		#if female selected
			
			if ($col[3] > $filter) {	

			open (OUT, ">>$outfile");
			print OUT "$line";
			$totalgenes++;			#log total genes

			}

		        else {

		        $filteredgenes++;			#log filtered genes
		        next

		        }
                }
			
	        elsif ($refsex eq "male") {             #if male selected

			if ($col[4] > $filter) {	

			open (OUT, ">>$outfile");
			print OUT "$line";
			$totalgenes++;			#log total genes

			}

		        else {

		        $filteredgenes++;			#log filtered genes
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

