#!/usr/bin/perl
#
# HOS_seq_extract.pl
#
# A script to report extract regions of sequences from Trace Archive files that might correspond to 
# tandem repeats with higher order structure (HOS)
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;
use FAlite;


###############################
#
#  Command-line options
#
################################

my $input; # specify input file
my $min_hos_level; # only look at HOS of a certain level
my $blast_db; # specify existing BLAST database
my $help; # print help

GetOptions ("input=s"         => \$input,
			"min_hos_level=i" => \$min_hos_level,
			"blast_db=s"      => \$blast_db,
		    "help"            => \$help
			);

###############################
#
# Set defaults
#
################################

# By default we are really most interested in looking for 3 levels of HOS, could change to 2 if 
# you want more results
$min_hos_level = 3  if (not defined $min_hos_level);


###############################
#
# Check command line options
#
################################

my $usage = "
usage: HOS_seq_extract.pl [options]
options:
  -input <file with HOS information (from HOS_finder.pl)>
  -min_hos_level (level of structure.) [$min_hos_level]
  -blast_db <BLAST database> - specify an existing BLAST database to use
  -help : print this help
";

die "Specify -input option to specify input file (hos.out.csv)\n" unless ($input);
die $usage if ($help);


###############################
#
# Variables
#
###############################

my $TMP = "/tmp";
my $output = "hos_repeats.fa";
unlink $output if -e $output;

##############################################################
#
# parse input file to get file names of sequences
#
##############################################################

# make a list of all files that we want to make our BLAST database from
my %files;

open(my $in, "<", $input) or die "Can't open $input\n";

# skip header
<$in>;

while(<$in>){
	chomp;
	my ($id, $hos_level, $file, $seq_id, $start, $end, $repeat_length, $copies) = split(/\t/);
	$files{$file} = 1 unless ($hos_level < $min_hos_level);
}
close($in);

my @files = keys %files;

my $BLASTDB = "HOS_BLASTDB";

# choose makeblastdb or formatdb (for older version of BLAST?) unless BLASTDB already exists
unless (defined($blast_db)){
	my $status = system("which makeblastdb > /dev/null;");
	if($status != 0){
		die "Can't find makeblastdb\n";
	} else {
		my $command = "makeblastdb -in \"@files\" -dbtype nucl -parse_seqids -out $BLASTDB";
		system("$command") && die "$command\n"; 
	}
}

# loop through input file again, this time to extract sequences
open($in, "<", $input) or die "Can't open $input\n";

# for each possible HOS repeat, we want to take note of the min and max coordinates 
my ($min, $max) = (0, 100000000);

# also want to keep track of when we move to each new HOS element in file
# i.e. track details of *previous* line in file
my $line_count = 0;
my %previous;

# skip header
<$in>;

while(<$in>){
	$line_count++;
	chomp;
	my ($id, $hos_level, $file, $seq_id, $start, $end, $repeat_length, $copies) = split(/\t/);

	if ($line_count == 1){
		$previous{id} = $id;
		$previous{seq_id} = $seq_id;
		$previous{start} = $start;
		$previous{end} = $end;
		$previous{repeat_lengths} = "$repeat_length,"; 
		$min = $start;
		$max = $end;
		next;
	}

	next if ($hos_level < $min_hos_level);
	
	# have we moved to a new HOS element?
	if ($id != $previous{id}){
		
		extract_seq($previous{seq_id}, $previous{id}, $min, $max, $previous{repeat_lengths});
			
		($min, $max) = (0, 10000000000);
		$previous{id} = $id;
		$previous{seq_id} = $seq_id;
		$previous{start} = $start;
		$previous{end} = $end;
		$previous{repeat_lengths} = "$repeat_length,"; 
		$min = $start;
		$max = $end;
	} else{
		# otherwise just modify $min and $max if necessary
		$min = $start if ($start < $min);
		$max = $end   if ($end   > $max);
		$previous{repeat_lengths} .= "$repeat_length,"; 
		
	}
}
close($in);

# once more for the last sequence in file
extract_seq($previous{seq_id}, $previous{id}, $min, $max, $previous{repeat_lengths});


exit(0);


#############################
#
#   S U B R O U T I N E S
#
#############################

sub extract_seq{
	my ($seq_id, $id, $start, $end, $repeat_lengths) = @_;
	$seq_id =~ s/\s+.*//;
	$repeat_lengths =~ s/,$//;
	
	my $command = "blastdbcmd -dbtype nucl -range $min-$max -db $BLASTDB -entry \"$seq_id\"";
#	print "$command\n";
	system("$command > $TMP/hos$$.fa") && die "Can't run $command\n";

	# now need to modify header of file and concatenate result to one main file
	# parse info from FASTA file
	open(my $in,  "<", "$TMP/hos$$.fa")     or die "Can't read from $TMP/hos$$.fa\n";
	open(my $out, ">", "$TMP/hos$$.fa.new") or die "Can't write to $TMP/hos$$.fa.new\n";

	my $fasta = new FAlite($in);

	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry){
	    my $seq = Keith::tidy_seq($entry->seq);
		my $header = $entry->def;
		$header =~ s/>//;
		my $new_header = ">HOS_ID=$id L=$repeat_lengths " . $header;
 		print $out "$new_header\n$seq\n";
	}
	close($in);
	close($out);
	
	system ("cat $TMP/hos$$.fa.new >> $output");
	unlink("$TMP/hos$$.fa");
	unlink("$TMP/hos$$.fa.new");    
}


