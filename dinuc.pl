#!/usr/bin/perl
#
# dinuc.pl 
#
# A script to count dinucleotide frequencies
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
	print "dinuc.pl - a program to calculate nucleotide and dinucleotide composition\n\n";
    print "1) By default, script will combine all sequences in input file before calculating nucleotide composition\n";
	print "2) Use -individual option to get statistics for individual sequences\n";
	print "4) -tsv mode prints one line of output per sequence in tab separated values, suitable for pasting into excel\n";
	print "5) Change precision of output with -precision option (default is 2 d.p.)\n";
    exit(0);
}


$precision = 2 if (!$precision);
			

# print help options if no arguments are specified							
if (!$ARGV[0]){
	system("$0 -help") && die "Could not run dinuc.pl -help\n";
	exit;
}



############################################################
# Read sequence file(s), add all sequence to $seq variable
############################################################

open(FILE,"$ARGV[0]") || die "Can't open $ARGV[0], please specify the FASTA sequence file you wish atcg.pl\nto analyse on the command line, e.g. dinuc.pl filename.\n\n";

my $fasta = new FAlite(\*FILE);

# need to store all sequences combined
my $combined_sequences;

# keep track of how many sequences there were
my $count = 0;

# print out header line, if we are using -tsv option
print "%AA\t%AC\t%AG\t%AT\t%CA\t%CC\t%CG\t%CT\t%GA\t%GC\t%GG\t%GT\t%TA\t%TC\t%TG\t%TT\n" if ($tsv);


# loop through each sequence in target file
while(my $entry = $fasta->nextEntry) {
	
	$count++;
		
	# are we calculating composition for individual sequences?
	if($individual){	
		# grab sequence counts and print output
		print_output(Keith::get_dinucleotides($entry->seq),$entry->def);
	}

	# add sequence to running total
	$combined_sequences .= ($entry->seq);
	
}
close(FILE);



# can now print output for combined sequence, unless using -tsv & -individual options
print_output(Keith::get_dinucleotides($combined_sequences),"All sequences combined") unless ($individual && $tsv);



################################################
#
# One main subroutine to handle printing output
#
################################################

sub print_output{
	my ($aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt,$o,$header) = @_;

	# first calculate total number of dinucleotides
	my $dinuc_count = $aa+$ac+$ag+$at+$ca+$cc+$cg+$ct+$ga+$gc+$gg+$gt+$ta+$tc+$tg+$tt+ $o;

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
	if ($dinuc_count == 0){

		# if using TSV output, just print zeroes as there is no sequence
		if($tsv){
			print "0\t" x 16,  "\n";
		}
		else{
			print "WARNING: no sequence detected\n\n\n\n";

		}
		return;
	}
	
	# calculate percentages
	my $aa_p = sprintf("%3.${precision}f",($aa / ($dinuc_count-$o))*100);
	my $ac_p = sprintf("%3.${precision}f",($ac / ($dinuc_count-$o))*100);
	my $ag_p = sprintf("%3.${precision}f",($ag / ($dinuc_count-$o))*100);
	my $at_p = sprintf("%3.${precision}f",($at / ($dinuc_count-$o))*100);
	my $ca_p = sprintf("%3.${precision}f",($ca / ($dinuc_count-$o))*100);
	my $cc_p = sprintf("%3.${precision}f",($cc / ($dinuc_count-$o))*100);
	my $cg_p = sprintf("%3.${precision}f",($cg / ($dinuc_count-$o))*100);
	my $ct_p = sprintf("%3.${precision}f",($ct / ($dinuc_count-$o))*100);
	my $ga_p = sprintf("%3.${precision}f",($ga / ($dinuc_count-$o))*100);
	my $gc_p = sprintf("%3.${precision}f",($gc / ($dinuc_count-$o))*100);
	my $gg_p = sprintf("%3.${precision}f",($gg / ($dinuc_count-$o))*100);
	my $gt_p = sprintf("%3.${precision}f",($gt / ($dinuc_count-$o))*100);
	my $ta_p = sprintf("%3.${precision}f",($ta / ($dinuc_count-$o))*100);
	my $tc_p = sprintf("%3.${precision}f",($tc / ($dinuc_count-$o))*100);
	my $tg_p = sprintf("%3.${precision}f",($tg / ($dinuc_count-$o))*100);
	my $tt_p = sprintf("%3.${precision}f",($tt / ($dinuc_count-$o))*100);
	
	# print concise output in TSV format, or more verbose details?
	if($tsv){
		print "$aa_p\t$ac_p\t$ag_p\t$at_p\t$ca_p\t$cc_p\t$cg_p\t$ct_p\t$ga_p\t$gc_p\t$gg_p\t$gt_p\t$ta_p\t$tc_p\t$tg_p\t$tt_p\n";
	}
	else{			
		print "Total number of dinucleotides: $dinuc_count\n\n";

		print "AA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tCG\tCT\tTA\tTC\tTG\tTT\tOther\n";
		print "$aa\t$ac\t$ag\t$at\t$ca\t$cc\t$cg\t$ct\t$ga\t$gc\t$gg\t$gt\t$ta\t$tc\t$tg\t$tt\t$o\n\n";

		print "%AA\t%AC\t%AG\t%AT\t%CA\t%CC\t%CG\t%CT\t%CA\t%CC\t%CG\t%CT\t%TA\t%TC\t%TG\t%TT\n";
		print "$aa_p\t$ac_p\t$ag_p\t$at_p\t$ca_p\t$cc_p\t$cg_p\t$ct_p\t$ga_p\t$gc_p\t$gg_p\t$gt_p\t$ta_p\t$tc_p\t$tg_p\t$tt_p\n\n";
		
	}
}


exit(0);

