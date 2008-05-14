#!/usr/bin/perl -w
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

use strict;
use FAlite;

open(IN, $ARGV[0]) or die "Can't open $ARGV[0]\n";
my $FA = new FAlite (\*IN);

open(OUT, ">$ARGV[0].rev") || die "Can't create output file\n";

while (my $entry = $FA->nextEntry) {
    my $seq = $entry->seq;
	
	my $def = $entry->def;
	my $revcom = reverse $seq;
	$revcom =~ tr/ACGTacgt/TGCAtgca/;
	my $tidied = &tidy_seq($revcom);
	print OUT "$def\n$tidied\n";
}

close(IN);
close(OUT);

sub tidy_seq{
#adds a new line character every 60 bases  
    my ($seq) = @_;
    $seq =~ s/[\s\n]//g;
    $seq =~ tr/a-z/A-Z/;
    
    my ($output_seq) = "";
    my (@seq2) = "";
    my ($start,$end);

    @seq2 = split(//,$seq);
    my $length = @seq2;
    my ($to_add) = int($length/60);

    $end = $start= 0;

    foreach (1..$to_add){
        $output_seq .= substr($seq,$start,60);
        $output_seq .= "\n";
        $start += 60;
    }
    $output_seq .= substr($seq,$start);
    return ($output_seq);
}

exit(0);
