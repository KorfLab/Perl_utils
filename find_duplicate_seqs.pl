#!/usr/bin/perl
#
# find_duplicate_seqs.pl
#
# A script to just report whether a file contains duplicate sequences
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;


# need a sequence file
if(!$ARGV[0]){
	die "Please specify the name of a valid FASTA sequence file\n";
}

    
# open seq file 
open(IN,"<$ARGV[0]") || die "Could not open $ARGV[0] file\n";

my $fasta = new FAlite(\*IN);


# two hashes, one to count sequences and one to store fasta header for each seq
my %seq2count;
my %seq2header;

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry) {

	# get and trim header
	my $header = $entry->def;
    $header =~ s/ .*//;
	
	# add info to hash
	my $seq = uc($entry->seq);
	
	$seq2count{$seq}++;
	$seq2header{$seq} = $header;
}

close(IN) || die "Couldn't close $ARGV[0]\n";


# now loop through each sequence to check the count
my $counter = 1;
foreach my $key (keys (%seq2count)){
	if($seq2count{$key}>1){
		
		# how many duplicates are there in total
		my $duplicates = $seq2count{$key} -1;
		
		
		print "$counter) The following entry has a sequence duplicated $duplicates time(s) in another entry\n";
		print "$seq2header{$key}\n";
		print "$key\n\n";
		$counter++;
	}
} 

exit(0);
