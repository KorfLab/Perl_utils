#!/usr/bin/perl
#
# trf_wrapper.pl
#
# A script to run the tandem repeat finder (trf) program and format the output
# acccording to our needs
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;


# options
my $match;
my $mismatch;
my $indel;
my $pmatch;
my $pindel;
my $min_score;
my $max_period;
my $min_copies; # minimum number of copies of repeat in read
my $min_length; # minimum length of repeat unit
my $min_repeat_fraction; # what percentage of the read length should be occupied by repeats?
my $high_repeat_threshold; # what fraction of read length will we use for our 'high repeat fraction' repeats

GetOptions ("match:i"   => \$match,
           "mismatch:i" => \$mismatch,
           "indel:i"    => \$indel,
           "pmatch:i"   => \$pmatch,
           "pindel:i"   => \$pindel,
           "min_score:i" => \$min_score,
           "period:i"   => \$max_period,
           "copies:f"   => \$min_copies,
           "length:i"   => \$min_length,
		   "min_repeat_fraction:f" => \$min_repeat_fraction,
		   "high_repeat_threshold:f" => \$high_repeat_threshold);

# set defaults if not specified on command line
# these are all defaults used (and suggested) by the trf program
$match = 1    if (!$match);
$mismatch = 1 if (!$mismatch);
$indel = 2    if (!$indel);
$pmatch = 80  if (!$pmatch);
$pindel = 5   if (!$pindel);
$min_score = 200  if (!$min_score); 
$max_period = 500 if (!$max_period);

# these are extra options that we can only implement through post-processing of raw trf output
$min_copies = 2 if (!$min_copies);
$min_length = 25  if (!$min_length);
$min_repeat_fraction = 0.5  if (!$min_repeat_fraction);
$high_repeat_threshold = 0.8 if (!$high_repeat_threshold);


# usage
die "
usage: trf.pl [options] <fasta file>
options:
  -match [$match]
  -mismatch [$mismatch]
  -indel [$indel]
  -pmatch (probability of match) [$pmatch]
  -pindel (probability of indel) [$pindel]
  -min_score [$min_score]
  -period (maximum period length) [$max_period]
  -copies (minimum copies) [$min_copies]
  -length (minimum repeat length) [$min_length]
  -min_repeat_fraction (minimum proportion of trace read that has to be repeat) [$min_repeat_fraction]
  -high_repeat_threshold (proportion of trace read that has to be repeat for HRF file) [$high_repeat_threshold]
" unless @ARGV == 1;

die "min_repeat_fraction must be lower than high_repeat_threshold\n" if ($min_repeat_fraction >= $high_repeat_threshold);

my ($file) = @ARGV;

my @param = ($match, $mismatch, $indel, $pmatch, $pindel, $min_score, $max_period);

# form name of output file from parameters
my $data_file  = $file . "." . join(".", @param) . ".dat";

# We may already have run trf on the same sequence file and just want to tweak the parameters for deciding which
# repeats to keep, so can make the script use an existing data file if one exists
if(-e $data_file){
	print "\nNOTE: A data file with the chosen parameters already exists ($data_file), will use this instead of re-running trf\n\n";
}
else{
	system("trf $file $match $mismatch $indel $pmatch $pindel $min_score $max_period -d -h > /dev/null") or die "Can't run trf\n";
#	system("trf $file $match $mismatch $indel $pmatch $pindel $min_score $max_period > /dev/null") or die "Can't run trf\n";	
}

# now need to capture lengths of all sequences in input file (frustratingly, this is not in the TRF output)
my %seq_to_length;
calculate_seq_lengths();

###############
# main loop
###############

# keep track of how many repeats are output
my $output_counter = 0;

# want to capture the (slightly processed) FASTA headers in the trf output
my $header;

# two output streams, one for repeats which make up a high repeat fraction (HRF) of the read
# and one for repeats which make up a low repeat fraction (LRF). These might contain repeat
# boundaries
open(OUT1,">$file.high.trf") or die "Can't open hrf output file\n";
open(OUT2,">$file.low.trf") or die "Can't open lrf output file\n";
open(DATA,$data_file) or die "Can't open data file\n";

my %errors;
my $seq_counter = 0;
my $repeat_counter = 0;

my $seq_repeat_counter; # to keep track of when a sequence contains more than one repeat
my @seq_repeat_data; # keep track of data on all repeats in a sequence to see if they are just multiples of each other

REPEAT: while(<DATA>){

	# skip blank likes
	next if (m/^$/);

	if (m/^Sequence: (.*)/) {
		$header = $1;
		$seq_counter++;

		# reset certain counters and data
		$seq_repeat_counter = 0;
		@seq_repeat_data = ();
		
		# and now move to next line in file
		next REPEAT;
	}
	
	# the main output that we are interested in will all be on one line which starts with various
	# numerical details of the repeat
	if(m/^\d+ \d+ \d+ \d+\.\d /){
		$repeat_counter++;
		
		# capture repeat data into various variables (most will be unsused)
		my ($start,$end,$period,$copies,$consensus,$matches,$indels,$score,$a,$c,$g,$t,$entropy,$repeat_seq) = split(/\s+/);
		
		my $repeat_length = length($repeat_seq);
		my $total_repeat_span = $end - $start + 1;
		my $repeat_fraction = $total_repeat_span / $seq_to_length{$header};

		my $tidied = Keith::tidy_seq($repeat_seq);
		

		###################################################################
		# now apply various criteria to see if we want to keep this repeat
		###################################################################
		
		# does the repeat occupy enough of the read that it is contained in?
		if ($repeat_fraction < $min_repeat_fraction){
			$errors{"min_repeat_fraction"}++;
			next REPEAT;
		}

		# is the repeat unit long enough? We certainly want to dicount tandem monomers, dimers etc.
		if ($repeat_length < $min_length){
			$errors{"repeat_length"}++;
			next REPEAT;
		}

		# are there enough copies of the repeat in the read?
		if ($copies < $min_copies){
			$errors{"min_copies"}++;
			next REPEAT;
		}

		# if we get this far then we will capture some of the repeat data in a hash 
		# this is to potentially compare to other repeats in the same sequence
		$seq_repeat_counter++;
		$seq_repeat_data[$seq_repeat_counter]{'coords'} = "$start $end";
		$seq_repeat_data[$seq_repeat_counter]{'consensus'} = "$consensus";
		$seq_repeat_data[$seq_repeat_counter]{'copies'} = "$copies";

		# Is this repeat just a multiple of a previously seen repeat in the same trace sequence?
		if ($seq_repeat_counter > 1){
			# loop through previous repeats
			for(my $i = 1; $i < @seq_repeat_data;$i++){

				if("$start $end" eq $seq_repeat_data[$i]{'coords'}){
					my $result = $consensus / $seq_repeat_data[$i]{'consensus'};
					my $processed_result = $result;
					$processed_result =~ s/\d+\.(\d+)/0\.$1/;
										
					# we are not expecting one repeat to be a perfect multiple of another repeat, so we'll allow 
					# some flexibility. 
					if(($processed_result < 0.15 || $processed_result > 0.85) && ($result > 1.5)){
#						print "$header\tRepeat $seq_repeat_counter: span = $start $end, consensus length = $consensus, copies = $copies\tTHIS IS A MULTIPLE ($result) OF:\n";
#						print "$header\tRepeat $i: span = $seq_repeat_data[$i]{'coords'}, consensus length = $seq_repeat_data[$i]{'consensus'}, copies = $seq_repeat_data[$i]{'copies'}\n";	
						$errors{"length_multiple"}++;
						next REPEAT;						
					}
				}
			}
		}
		
		# if we are here then we are keeping the repeat and printing it out
		$output_counter++;

		my $formatted = sprintf("%.0f",$repeat_fraction * 100);
		# which output file will the repeats go in?
		if($repeat_fraction > $high_repeat_threshold){
			print OUT1 ">tandem-$output_counter N=$copies L=$repeat_length F=$formatted% P=$header\n";
			print OUT1 "$tidied\n";			
		}
		else{
			print OUT2 ">tandem-$output_counter N=$copies L=$repeat_length F=$formatted% P=$header\n";
			print OUT2 "$tidied\n";
		}
	}
}
close(OUT1);
close(OUT2);
close(DATA);


###############################
#
# Print summary statistics
#
################################

print "Processed $seq_counter sequences containing $repeat_counter repeats";
print ", $output_counter of which matched all criteria\n\n";

my $total_rejections = $errors{'repeat_length'} + $errors{'min_copies'} + $errors{'min_repeat_fraction'} + $errors{'length_multiple'};
print "$total_rejections repeats rejected for failing criteria:\n";
print "\t$errors{'min_repeat_fraction'} repeats rejected for making up too low a fraction of the total sequence (<$min_repeat_fraction)\n";
print "\t$errors{'repeat_length'} repeats rejected for being too short (<$min_length nt)\n";
print "\t$errors{'min_copies'} repeats rejected for having too few copies (<$min_copies)\n";
print "\t$errors{'length_multiple'} repeats rejected for being a multiple of a shorter repeat\n\n";


exit(0);


###############################
#
#   S U B R O U T I N E S
#
################################

sub calculate_seq_lengths{
	open (FILE,"<$file") || die "Can't open file $file\n";
	my $fasta = new FAlite(\*FILE);

	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry) {

		# extract seq length & header and add to hash
	    my $length = length($entry->seq);		
		my $header = $entry->def;
		$header =~ s/>//;
		$seq_to_length{$header} = $length;
		
	}
	close(FILE);
}
