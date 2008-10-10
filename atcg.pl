#!/usr/bin/perl
#
# atcg.pl 
#
# A script to count GC content plus mononucleotide frequencies
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use FAlite;
use Keith;

my $help;       # print help
my $tsv;	    # print main part of output suitable for parsing by a program or excel (tab separated values)
my $individual; # print out individual statistics for each sequence (if there are individualple)
my $precision;  # how many decimal places to use for output

GetOptions ("help"        => \$help,
			"tsv"	      => \$tsv,
			"individual"  => \$individual,
			"precision=i" => \$precision);


if($help){
	print "atcg.pl - a program to calculate nucleotide and dinucleotide composition\n\n";
    print "1) By default, script will combine all sequences in input file before calculating nucleotide composition\n";
	print "2) Use -individual option to get statistics for individual sequences\n";
    print "3) GC content is calculated both for the entire sequence, and also for a masked copy of the sequence\n";
    print "   (calculated after excluding any Ns or non-ACGT characters)\n";
	print "4) -tsv mode prints one line of output per sequence in tab separated values, suitable for pasting into excel\n";
	print "5) Change precision of output with -precision option (default is 2 d.p.)\n";
    exit(0);
}


$precision = 2 if (!$precision);
			

# print help options if no arguments are specified							
if (!$ARGV[0]){
	system("$0 -help") && die "Could not run atcg.pl -help\n";
	exit;
}



############################################################
# Read sequence file(s), add all sequence to $seq variable
############################################################

open(FILE,"$ARGV[0]") || die "Can't open $ARGV[0], please specify the FASTA sequence file you wish atcg.pl\nto analyse on the command line, e.g. atcg.pl filename.\n\n";

my $fasta = new FAlite(\*FILE);

# need to store all sequences combined
my $combined_sequences;

# keep track of how many sequences there were
my $count = 0;

# print out header line, if we are using -com option
print "%A\t%C\t%G\t%T\t%N\t%Other\tLength\t%GC\t%GC (masked)\n" if ($tsv);


# loop through each sequence in target file
while(my $entry = $fasta->nextEntry) {
	
	$count++;
		
	# are we calculating composition for individual sequences?
	if($individual){	
		# grab sequence counts and print output
		print_output(Keith::get_mononucleotides($entry->seq),$entry->def);
	}

	# add sequence to running total
	$combined_sequences .= ($entry->seq);
	
}
close(FILE);



# can now print output for combined sequence, unless using -tsv & -individual options
print_output(Keith::get_mononucleotides($combined_sequences),"All sequences combined") unless ($individual && $tsv);



################################################
#
# One main subroutine to handle printing output
#
################################################

sub print_output{
	my ($a,$c,$g,$t,$n,$o,$header) = @_;

	# first calculate length
	my $total_length = $a + $c + $g + $t + $n + $o;

	# print out header if we are individual mode, else print total sequence count
	# but don't print this out with -com option selected
	if($individual && !$tsv){
		print '-' x 60, "\n";
		print "$header\n";
		print '-' x 60, "\n";
	}
	elsif(!$tsv){
		print "\nProcessed $count sequences\n\n";
	}
	
	# if there was no sequence (for some reason), we should print an error and stop
	if ($total_length == 0){

		# if using TSV output, just print zeroes as there is no sequence
		if($tsv){
			print "0\t" x 9,  "\n";
		}
		else{
			print "WARNING: no sequence detected\n\n\n\n";

		}
		return;
	}
	

	# Calculate final statistics
	my $gc_percent  = sprintf("%3.${precision}f",(($g + $c)/$total_length)*100);

	# masked length ignores Ns (or other non-AGCT characters) in sequence
	my $masked_length     = $a + $c + $g + $t;
	my $masked_gc_percent = sprintf("%3.${precision}f",(($g + $c)/$masked_length)*100);

	my $a_percent         = sprintf("%3.${precision}f",($a / $total_length)*100);
	my $c_percent         = sprintf("%3.${precision}f",($c / $total_length)*100);
	my $g_percent         = sprintf("%3.${precision}f",($g / $total_length)*100);
	my $t_percent         = sprintf("%3.${precision}f",($t / $total_length)*100);
	my $n_percent         = sprintf("%3.${precision}f",($n / $total_length)*100);
	my $o_percent         = sprintf("%3.${precision}f",($o / $total_length)*100);
	
	# print concise output in TSV format, or more verbose details?
	if($tsv){
		print "$a_percent\t$c_percent\t$g_percent\t$t_percent\t$n_percent\t$o_percent\t$total_length\t$gc_percent\t$masked_gc_percent\n";
	}
	else{			
		print "Length\tGC%\n";
		print "$total_length\t$gc_percent\n\n";

		print "Masked length\tMasked GC%\n";
		print "$masked_length\t$masked_gc_percent\n\n";

		print "A\tC\tG\tT\tN\tOther\n";
		print "$a\t$c\t$g\t$t\t$n\t$o\n\n";

		print "A%\tC%\tG%\tT%\tN%\tOther%\n";
		printf "%3.${precision}f\t%3.${precision}f\t%3.${precision}f\t%3.${precision}f\t%3.${precision}f\t%3.${precision}f\n\n\n\n", 
		$a_percent, $c_percent, $g_percent, $t_percent, $n_percent, $o_percent;
		
	}
}


exit(0);

