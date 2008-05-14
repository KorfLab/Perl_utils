#!/usr/bin/perl -w
#
# tidy.pl
#
# simple script to tidy a FASTA file to 60 characters per line
# and make sequence all the same case
#
# by Keith Bradnam
#
# created 11th November 2005
#
# Last updated by: $Author$
# Last updated on: $Date$    
#
#######################################################################################

use strict;

unless ($ARGV[0]){
    print"please enter the input file name: ";
    chop ($ARGV[0]=<STDIN>);
}
until (open (INFILE,$ARGV[0])){
    print "sorry, $ARGV[0] was not found, ";
    print"please re-enter the file name: ";
    chop ($ARGV[0]=<STDIN>);
}

unless ($ARGV[1]){
    print"please enter the output file name: ";
    chop ($ARGV[1]=<STDIN>);
}


open (OUT,">$ARGV[1]")|| die "$!";


##################################################

my ($flag,$seq);

SEQ:while(<INFILE>){
    if(/^>/){
		if($flag){
	   		print OUT &tidy_seq($seq),"\n";
		}
		print OUT;
		$seq = "";
		$flag = 1;
		next SEQ;
    }
    $seq .= $_;
}

if($flag){
   	print OUT &tidy_seq($seq),"\n";
}


print "\n";				
###########################################################
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

######################################
