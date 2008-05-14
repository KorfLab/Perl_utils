#!/usr/bin/perl
#
# seqs2size.pl
#
# a quick script to calculate the size of sequences in a
# fasta file
#
# by Keith Bradnam, 14th September 2005
##############################################

use strict;
use warnings;
use List::Util qw(sum);
use Getopt::Long;
use FAlite;

# command line options
# just one: -summary mode prints summary statistics of all sequences in a file
my $summary; 
GetOptions ("summary" => \$summary);


# array to store all sequence lengths
my @sizes;

open (FILE,"<$ARGV[0]") || die "Can't open file $ARGV[0]\n";
my $fasta = new FAlite(\*FILE);

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry) {
	
	my $seq = $entry->seq;
    my $length = length($seq);
	
	# are we counting stats or just printing sequence lengths?
	if ($summary){
		push(@sizes,$length);
	}
	else{
		print "$length\n";	
	}
}
close(FILE);

# can quit unless we are in summary mode
exit(0) unless $summary; 

# calculate and print statistics
my $n = scalar(@sizes);
my $mean = sprintf("%.2f",sum(@sizes)/$n);

# no point continuing if sample size is 1
die "Can't calculate meaningful statistics on 1 sequence\n" if ($n == 1);
    
my $stdev = sprintf("%.2f",sqrt(sum(map {($_ - $mean) ** 2} @sizes) / ($n-1)));

# standard error
 my $se = sprintf("%.2f",$stdev / sqrt($n));

# 95% confidence limits of the mean
my $ci = $se * 1.96;

# warn if sample size is too low to use 1.96 for calculating confidence intervals
print "Sample size ($n) is low, calculation of 95% confidence limits will be underestimates\n\n" if ($n < 100);

# print final stats.
print "\n";
print "N\tMean\tStd_dev\tStd_err\t95%_CI\n";
print "$n\t$mean\t$stdev\t$se\t$ci\n\n";


exit(0);