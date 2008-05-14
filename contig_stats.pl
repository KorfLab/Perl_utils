#!/usr/bin/perl -w
#
# contig_stats.pl
#
# a quick script to calculate the mean and N50 sizes for a file of
# contig or scaffold sequences
#
# by Keith Bradnam
#
# Last updated by: $Author$
# Last updated on: $Date$
##############################################

use strict;

my $flag;
my $seq;
my $length;			# for individual contigs
my @contig_sizes;   # lengths stored in this array 
my $total_length; 	# sum of all contig lengths
my $seq_count;		# count number of contigs
my $n;				# keep track of N's in sequence

SEQ: while(<>){
	chomp;
  	if(/^>/){
		$seq_count++;
		if($flag){
	     	$length = length($seq);
		 	push(@contig_sizes, $length);
		 	$total_length += $length;
		 	# now calculate N's in sequence (two passes to capture both lower and upper case)
		 	$n +=  $seq =~ tr/N/N/;
		 	$n += $seq =~ tr/n/n/;		 	
	 	}
		$seq = "";
		$flag = 1;
    	next SEQ;
    	}
   $seq .= $_;
}

# get last sequence in file
$length = length($seq);
push(@contig_sizes, $length);
$total_length += $length;

# now calculate N's in sequence (two passes to capture both lower and upper case)
$n += $seq =~ tr/N/N/;
$n += $seq =~ tr/n/n/;		 	



# calculate N's as a percentage
my $percent_N = sprintf("%.1f",($n/$total_length)*100);

print "\n$total_length - total assembly length (bp)\n";
print "$n - N's (${percent_N}%)\n";
print "$seq_count - number of contigs\n";
my $mean_size = sprintf("%.1f",$total_length/$seq_count);
print "$mean_size - average contig size (bp)\n";


# Now need to loop through (sorted) list of contig sizes
# N50 size is the length 'x' such that 50% of the sequence is contained in contigs 
# of length 'x' or greater
my $running_total;
foreach my $contig(reverse sort{$a <=> $b;} (@contig_sizes)){
	$running_total+=$contig;
	if($running_total/$total_length >=.5){
		print "$contig - N50 length (bp)\n\n";
		last;
	}    
}

exit(0);