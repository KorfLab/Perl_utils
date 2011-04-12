#!/usr/bin/perl
#
# price2trf.pl
#
# A script to run take assembled contigs from the PRICE assembler and produce a set
# of tandem repeats from them (using trf). Based on trf_wrapper.pl
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;


###############################
#
#  Command-line options
#
################################

my $match;
my $mismatch;
my $indel;
my $pmatch;
my $pindel;
my $min_score;
my $max_period;
my $min_copies; # minimum number of copies of repeat in read
my $min_length; # minimum length of repeat unit
my $stats;      # show summary statistics fir repeats that are skipped/ignored
my $help; # print help

GetOptions ("match:i"     => \$match,
            "mismatch:i"  => \$mismatch,
            "indel:i"     => \$indel,
            "pmatch:i"    => \$pmatch,
			"pindel:i"    => \$pindel,
            "min_score:i" => \$min_score,
            "period:i"    => \$max_period,
            "copies:f"    => \$min_copies,
            "length:i"    => \$min_length,
			"stats"       => \$stats,
		    "help"        => \$help		
			);


###############################
#
# Set defaults
#
################################

# set defaults if not specified on command line
# these are all defaults used (and suggested) by the trf program
$match = 1    if (!$match);
$mismatch = 1 if (!$mismatch);
$indel = 2    if (!$indel);
$pmatch = 80  if (!$pmatch);
$pindel = 5   if (!$pindel);
$min_score = 200  if (!$min_score); 
$max_period = 750 if (!$max_period);

# these are extra options that we can only implement through post-processing of raw trf output
$min_copies = 2   if (!$min_copies);
$min_length = 5  if (!$min_length);


###############################
#
# Check command line options
#
################################

# usage
my $usage = "
usage: trf_wrapper.pl [options] <input FASTA file>
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
  -stats : show details of which tandem repeats were skipped or ignored
  -help : print this help
";

die $usage if ($help);
die $usage unless (@ARGV == 1);
my ($file) = @ARGV;

###############################
#
# Variables
#
###############################

my @param = ($match, $mismatch, $indel, $pmatch, $pindel, $min_score, $max_period);

my %seq_to_length; # need to know lengths of all sequences, key to hash is the fasta header
my %errors; # Keep track of various errors
$errors{'repeat_length'} = $errors{'min_copies'} = $errors{'min_repeat_fraction'} = $errors{'length_multiple'} = 0;

my $seq_counter = 0;      # how many input sequences do we process
my $repeat_counter = 0;   # number of tandem repeats predicted by trf
my $output_counter = 0;   # number of tandem repeats that pass all filters 
my $total_nt = 0;         # keep track of the total size of sequence in the input files
my $repeat_cutoff = 0.25; # what fraction of a contig has to be made up on tandem repeats in order for us to keep it
my $excluded_trf_nt = 0;  #


###############################
#
# Loop over input file
#
###############################

# form name of output file from parameters
my $data_file  = $file . "." . join(".", @param) . ".dat";

# We may already have run trf on the same sequence file and just want to tweak the parameters for deciding which
# repeats to keep, so can make the script use an existing data file if one exists
if(-e $data_file){
	warn "NOTE: A data file with the chosen parameters already exists ($data_file), will use this instead of re-running trf\n\n";
}
else{
	system("trf $file $match $mismatch $indel $pmatch $pindel $min_score $max_period -d -h > /dev/null") or die "Can't run trf\n";
}		

# now process trf output file
process_trf_output($file, $data_file);



# Statistics on what was rejected/outside of thresholds
if($stats){
	my $total_rejections = $errors{'repeat_length'} + $errors{'min_copies'} + $errors{'min_repeat_fraction'} + $errors{'length_multiple'};
	warn "\nProcessed $seq_counter sequences that contained $repeat_counter repeats, $output_counter of which matched all criteria\n";
	warn "$total_rejections repeats rejected for failing criteria:\n";
	warn "\t$errors{'min_repeat_fraction'} repeats rejected for making up too low a fraction of the total sequence (<$repeat_cutoff)\n";
	warn "\t$errors{'repeat_length'} repeats rejected for being too short (<$min_length nt)\n";
	warn "\t$errors{'min_copies'} repeats rejected for having too few copies (<$min_copies)\n";
	warn "\t$errors{'length_multiple'} repeats rejected for being a multiple of a shorter repeat\n\n";

	my $excluded_percent  = sprintf("%.2f",($excluded_trf_nt / $total_nt)*100);
	warn "Total amount of tandem repeats in excluded repeat fraction (<$repeat_cutoff): ", int($excluded_trf_nt), " nt ($excluded_percent%)\n\n";	
}

exit(0);



###############################
#
#
#   S U B R O U T I N E S
#
#
################################


sub process_trf_output{

	my ($file,$data) = @_;
	
	# want to capture the (slightly processed) FASTA headers in the trf output
	my $header;

	
	# now need to capture lengths of all sequences in input file (frustratingly, this is not in the TRF output)
	calculate_seq_lengths($file);
	
	open(DATA,"<$data") or die "Can't open data file\n";

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
			if ($repeat_fraction < $repeat_cutoff){
				$errors{"min_repeat_fraction"}++;
				$excluded_trf_nt += ($copies * $repeat_length);
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
			$seq_repeat_data[$seq_repeat_counter]{'coords'}    = "$start $end";
			$seq_repeat_data[$seq_repeat_counter]{'consensus'} = "$consensus";
			$seq_repeat_data[$seq_repeat_counter]{'copies'}    = "$copies";

			# Is this repeat just a multiple of a previously seen repeat in the same trace sequence?
			if ($seq_repeat_counter > 1){
				# loop through previous repeats
				for(my $i = 1; $i < @seq_repeat_data; $i++){

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
			my $mass = $repeat_length * $copies; 
			print ">tandem-$output_counter L=$repeat_length N=$copies M=$mass F=$formatted% P=$header\n";
			print "$tidied\n";			
		}
	}
	close(DATA);
}



sub calculate_seq_lengths{
	my ($file) = @_;
	open (FILE,"<$file") || die "Can't open file $file\n";
	my $fasta = new FAlite(\*FILE);

	while(my $entry = $fasta->nextEntry) {
		# extract seq length & header and add to hash and to grand total
	    my $length = length($entry->seq);
		$total_nt += $length;
		my $header = $entry->def;
		$header =~ s/>//;
		$seq_to_length{$header} = $length;		
	}
	close(FILE);
}
