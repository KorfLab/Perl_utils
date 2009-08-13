#!/usr/bin/perl
#
# entropy.pl 
#
# A script to count the entropy of a set of DNA sequences
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use FAlite;
use Keith;

# print help options if no arguments are specified							
if (@ARGV != 1){
	die "Usage: entropy.pl <filename>\n";
}



############################################################
# Read sequence file(s), add all sequence to $seq variable
############################################################

open(FILE,"$ARGV[0]") || die "Can't open $ARGV[0]\n";

my $fasta = new FAlite(\*FILE);

# need to store all sequences combined
my $combined_sequences;

# keep track of how many sequences there were
my $count = 0;

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry) {
	
	$count++;

	# add sequence to running total
	$combined_sequences .= ($entry->seq);
	
}
close(FILE);


my ($a,$c,$g,$t,$n,$o) = Keith::get_mononucleotides($combined_sequences);

# first calculate sum of just a, c, g, and t (there may be N's and other nucleotides)
my $total_length = $a + $c + $g + $t;

# if there was no sequence (for some reason), we should print an error and stop
if ($total_length == 0){
	die "WARNING: no sequence detected\n";
}

my $precision = 4;

my $a_freq = sprintf("%3.${precision}f",($a / $total_length));
my $c_freq = sprintf("%3.${precision}f",($c / $total_length));
my $g_freq = sprintf("%3.${precision}f",($g / $total_length));
my $t_freq = sprintf("%3.${precision}f",($t / $total_length));

my $entropy = 0;
$entropy += ($a_freq * (log($a_freq)/log(2)));
$entropy += ($c_freq * (log($c_freq)/log(2)));
$entropy += ($g_freq * (log($g_freq)/log(2)));
$entropy += ($t_freq * (log($t_freq)/log(2)));

$entropy = sprintf("%3.${precision}f",-($entropy));

print "Entropy of $count sequences is $entropy\n";

exit(0);

