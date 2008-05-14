#!/usr/bin/perl
#
# atcg.pl 
#
# A script to count GC content plus mono- and dinucleotide frequencies
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use FAlite;

my $help;    # print help
my $dinuc;   # calculate dinucleotide frequencies as well as mononucleotide frequencies
my $com;	 # print main part of output suitable for parsing by a program or excel (computer readable)

GetOptions ("dinuc"   => \$dinuc,
	    	"help"    => \$help,
			"com"	  => \$com);


if($help){
	print "atcg.pl - a program to calculate nucleotide and dinucleotide composition\n\n";
    print "1) atcg.pl can analyse multiple sequences, but will concatenate the results together.\n";
    print "2) use -dinuc option to also calculate dinucleotide frequencies\n";
    print "3) GC content is calculated irrespective of Ns in input sequence, Masked GC% is\n";
    print "   calculated using length of sequence MINUS any N characters\n";
	print "4) -com mode prints one line output of most of the data, suitablef or pasting into excel\n\n";
    exit(0);
}


my %mono = ('A' => 0,'C' => 0,'G' => 0,'T' => 0,'N' => 0);

# hashes for dinucleotide frequencies and percentages
my %di;
my %di_percent;

# count sequences in input file
my $n = 0;
							
############################################################
# Read sequence file(s), add all sequence to $seq variable
############################################################

open(FILE,"$ARGV[0]") || die "Can't open $ARGV[0], please specify the FASTA sequence file you wish atcg.pl\nto analyse on the command line, e.g. atcg.pl filename.\n\n";


my $fasta = new FAlite(\*FILE);

# loop through each sequence in target file
while(my $entry = $fasta->nextEntry) {	
	$n++;
	my $seq = uc($entry->seq);

	# count mononucleotide frequencies
	$mono{"A"} += ($seq =~ tr/A/A/); 
	$mono{"C"} += ($seq =~ tr/C/C/); 
	$mono{"G"} += ($seq =~ tr/G/G/); 
	$mono{"T"} += ($seq =~ tr/T/T/); 
	$mono{"N"} += ($seq =~ tr/N/N/); 

	# do we also want to calculate dinucleotides? Skip if we don't
	next unless ($dinuc || $com);
	
	# to count dinucleotides, loop through sequence, take 2 bp and increment the hash counter
	foreach my $i (0..length($seq)){
	    my $tmp = substr($seq,$i,2);		
		$di{$tmp}++;
	}
}
close(FILE);

# Calculate final statistics
# need to determine number of Ns to work out masked length and GC content of masked sequence
my $total_length = $mono{'A'} + $mono{'C'} + $mono{'G'} + $mono{'T'} + $mono{'N'};
my $masked_length     = $total_length - $mono{'N'};
my $total_GC_percent  = sprintf("%3.2f",(($mono{'G'}+$mono{'C'})/$total_length)*100);
my $masked_GC_percent = sprintf("%3.2f",(($mono{'G'}+$mono{'C'})/$masked_length)*100);
my $A_percent         = ($mono{'A'}/$total_length)*100;
my $T_percent         = ($mono{'T'}/$total_length)*100;
my $C_percent         = ($mono{'C'}/$total_length)*100;
my $G_percent         = ($mono{'G'}/$total_length)*100;
my $N_percent         = ($mono{'N'}/$total_length)*100;

# print human or computer readable output?
unless ($com){
	print "\nProcessed $n sequences\n\n";
	print "Length\tGC%\tMasked length\tMasked GC%:\n";
	print "$total_length\t$total_GC_percent\t$masked_length\t$masked_GC_percent\n\n";

	print "A\tC\tG\tT\tN\n";
	print "$mono{'A'}\t$mono{'C'}\t$mono{'G'}\t$mono{'T'}\t$mono{'N'}\n\n";

	print "A%\tC%\tG%\tT%\tN%\n";
	printf "%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n\n", $A_percent, $C_percent, $G_percent, $T_percent, $N_percent;
}

exit(0) unless ($dinuc || $com);



######################################
# calculate dinucleotide percentages
######################################

				

my $dinuc_total = $di{'AA'}+$di{'AC'}+$di{'AG'}+$di{'AT'} + 
				  $di{'CA'}+$di{'CC'}+$di{'CG'}+$di{'CT'} +
				  $di{'GA'}+$di{'GC'}+$di{'GG'}+$di{'GT'} +
				  $di{'TA'}+$di{'TC'}+$di{'TG'}+$di{'TT'};

# will make four lines of output text
my ($out1,$out2,$out3,$out4);				
foreach my $i qw(AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT){
	
	# set to 0 if no values exist for that dinucleotide
	($di{$i} = 0) if (!defined($di{$i}));
	
	$di_percent{$i} = (($di{$i}/$dinuc_total)*100);
	$out1 .= "$i\t";
	$out2 .= "$di{$i}\t";
	$out3 .= "$i%\t";
	$out4 .= sprintf("%3.2f\t",$di_percent{$i});
}

if(!$com){
	print "$out1\n$out2\n\n$out3\n$out4\n\n";
}
else{
	print "\nN\tLength\tGC%\tA%\tC%\tG%\tT%\tN%\t$out3\n";
	print "$n\t$total_length\t$total_GC_percent\t";
	printf "%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t", $A_percent, $C_percent, $G_percent, $T_percent, $N_percent;
	print "$out4\n\n";

}

exit(0);

