#!/usr/bin/perl
#
# revcomp.pl
#
# simple script to reverse complement a FASTA sequence file
# also tidies it to make it 60 characters per line
#
# by Keith Bradnam
#
# created 9th July 2007
#
# Last updated by: $Author$
# Last updated on: $Date$    
#
#######################################################################################

use strict; use warnings;
use Keith;
use FAlite;

die "Usage: revcomp.pl <input_FASTA_file>\n" unless (@ARGV == 1);

open(IN, $ARGV[0]) or die "Can't open $ARGV[0]\n";
my $FA = new FAlite (\*IN);

open(OUT, ">$ARGV[0].rev") || die "Can't create output file\n";

while (my $entry = $FA->nextEntry) {
    my $seq = $entry->seq;
	
	my $def = $entry->def;
	my $revcom = reverse $seq;
	$revcom =~ tr/ACGTacgt/TGCAtgca/;
	my $tidied = Keith::tidy_seq($revcom);
	print OUT "$def\n$tidied\n";
}

close(IN);
close(OUT);


exit(0);
