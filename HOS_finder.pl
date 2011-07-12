#!/usr/bin/perl
#
# HOS_finder.pl
#
# A script to report on possible Higher Order Structure of tandem repeats present in Trace Archive reads
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

my $dir; # use current directory for data files
my $min_copies; # minimum number of copies of repeat in read
my $min_length; # minimum length of repeat unit
my $offset; # what offset between overlapping repeats in allowed 
my $help; # print help

GetOptions ("dir=s"    => \$dir,
            "copies:f" => \$min_copies,
            "length:i" => \$min_length,
			"offset:i" => \$offset, 
		    "help"     => \$help
			);

###############################
#
# Set defaults
#
################################

# these are extra options that we can only implement through post-processing of raw trf output
$min_copies = 2  if (not defined $min_copies);
$min_length = 25 if (not defined $min_length);
$offset = 5      if (not defined $offset);

###############################
#
# Check command line options
#
################################

my $usage = "
usage: HOS_finder.pl [options]
options:
  -dir <directory name containing fasta files>
  -copies (minimum copies) [$min_copies]
  -length (minimum repeat length) [$min_length]
  -offset (offset allowed between overlapping repeats) [$offset]
  -help : print this help
";

die "Specify -dir option to process all *.dat files in specified directory\n" unless ($dir);
die $usage if ($help);


###############################
#
# Variables
#
###############################


# hash to store sequence lengths
my %seq_to_length;

# one hash to store details of repeats with higher order structure (at least 3 different repeat lengths)
my %hos;


###############################
#
# Loop through input file(s)
#
###############################

print STDERR "\n# $0 started at ", `date`, "\n";

foreach my $file (glob("$dir/*processed_traces*fa*dat")){
	print STDERR "Processing $file\n";
	
	# skip zero byte files
	next unless (-s $file);
			
	# want to capture the (slightly processed) FASTA headers in the trf output
	my $header;
	
	open(my $data,"<", "$file") or die "Can't open $file\n";

	my $seq_repeat_counter; # to keep track of when a sequence contains more than one repeat
	my @seq_repeat_data;    # to keep track of data on all repeats in a sequence to see if they are just multiples of each other

	# need a generic counter for each sequence that might contain tandem repeats
	# will form part of a unique ID later on
	my $seq_count = 0;
	
	REPEAT: while(<$data>){

		# skip blank likes
		next if (m/^$/);

		# Need to know when we get to each new sequence in file (1 sequence might contain multiple tandem repeats)
		if (m/^Sequence: (.*)/) {
			$header = $1;
			
			$seq_count++;
			
			# reset certain counters and data
			$seq_repeat_counter = 0;
			@seq_repeat_data = ();

			# and now move to next line in file
			next REPEAT;
		}

		# the main output that we are interested in will all be on one line which starts with various
		# numerical details of the repeat
		if(m/^\d+ \d+ \d+ \d+\.\d /){

			# capture repeat data into various variables (most will be unsused)
			my ($start,$end,$period,$copies,$consensus,$matches,$indels,$score,$a,$c,$g,$t,$entropy,$repeat_seq) = split(/\s+/);

			###################################################################
			# now apply various criteria to see if we want to keep this repeat
			###################################################################

			# are there enough copies of the repeat in the read? N.B. you can have less than 2 copies present in TRF
			next REPEAT unless ($copies >= $min_copies);

			# is the repeat length at least $min_length (we want to avoid comparing things like 8 nt vs 4 nt repeats)
			next REPEAT unless ($consensus >= $min_length);

			# if we get this far then we will capture some of the repeat data in a hash 
			# this is to potentially compare to other repeats in the same sequence
			$seq_repeat_counter++;
			$seq_repeat_data[$seq_repeat_counter]{start}     = "$start";
			$seq_repeat_data[$seq_repeat_counter]{end}       = "$end";
			$seq_repeat_data[$seq_repeat_counter]{consensus} = "$consensus";
			$seq_repeat_data[$seq_repeat_counter]{copies}    = "$copies";
			$seq_repeat_data[$seq_repeat_counter]{file}      = "$file";
			
			# add repeat info to %HOS hash, but only if it is first repeat in sequence
			my $key = "$header,$seq_count";
			my $info = "$start,$end,$consensus,$copies,$file,1";
			push(@{$hos{$key}}, $info) if ($seq_repeat_counter == 1);
						
			# no point going any further until we have seen at least 2 repeats within the current sequence
			next REPEAT unless ($seq_repeat_counter > 1);
			
			# loop through previous repeats. Want to see if one is a multiple of another
			for (my $i = 1; $i < $seq_repeat_counter; $i++){

				# want to see if current start/end coordinates are identical or within $offset nt of previous repeat
				next unless (abs($start - $seq_repeat_data[$i]{start}) <= $offset);
				next unless (abs($end   - $seq_repeat_data[$i]{end})   <= $offset);

				# $result is the ratio of the longest repeat to the shorter one, need longer repeat to be at least 1.85 x length of shorter repeat
				# how we calculate this depends on whether current repeat length is longer than previous one or not
				my $result;
				
				if ($consensus < $seq_repeat_data[$i]{consensus}){
					$result = $seq_repeat_data[$i]{'consensus'}/ $consensus;
				} else{
					$result = $consensus / $seq_repeat_data[$i]{'consensus'};
				}
				my $processed_result = $result;
				$processed_result =~ s/\d+\.(\d+)/0\.$1/;

				next unless (($processed_result < 0.15 or $processed_result > 0.85) && ($processed_result > 1.5));
				
				# we should now be looking at a multiple repeat so add it to hash
				my $info = "$start,$end,$consensus,$copies,$file,$result";
				push(@{$hos{$key}}, $info);
								
				# if we have seen a match, don't need to look any further
				next REPEAT;
			}
		}
	}
	close($data);
}


# want a hash to keep track of how many 2nd-level, 3rd-level etc. HOS are preesnt
my %hos_count;

# now loop through %hos structure and write outptut to file
open(my $out, ">", "hos.out.csv") or die "Can't write to hos.out.csv\n";

print $out "HOS_ID\tHOS_LEVEL\tFILE\tSEQ_ID\tSTART\tEND\tREPEAT_LENGTH\tCOPIES\n";
my $counter = 1;
foreach my $key (keys %hos){
	my $matches = @{$hos{$key}};
	$hos_count{$matches}++;
	my ($seq_id, $seq_count) = split(/,/,$key);
	if ($matches > 2){
		foreach my $match (@{$hos{$key}}){
			my ($start, $end, $repeat_length, $copies, $file, $n) = split(/,/, $match);
			$file =~ s#^\.\/##;
			$file =~ s/(.*\.\d{3}\.fa).*/$1/;
			print $out "$counter\t$matches\t$file\t$seq_id\t$start\t$end\t$repeat_length\t$copies\n";
		}
	$counter++;
	}
}
close($out);
print "\n";
foreach my $key (sort {$a <=> $b} keys %hos_count){
	print "HOS level $key -> $hos_count{$key} occurences found\n";
}
print STDERR "\n# $0 finished at ", `date`, "\n";

exit(0);


