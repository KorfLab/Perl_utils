#!/usr/bin/perl
#
# ftp_trace_reads.pl
#
# A script to download trace reads of selected species from NCBI trace archive
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use Net::FTP;
use Cwd;
use List::Util 'shuffle';


###############################
#
#  Set some environment variables
#
################################

$SIG{'INT'} = 'INT_handler'; # capture signal interrupts
$ENV{FTP_PASSIVE} = 1;       # Need this or else FTP won't work in (the default) active mode!


###############################
#
#  Command-line options
#
################################

my $max_files;        # maximum number of files to download for any species
my $debug;            # whether to turn on debugging in FTP module
my $timeout;          # set timeout value for Net::FTP module
my $sleep;            # how long to sleep for before retrying ftp download
my $max_attempts;     # how many attempts to download one file before giving up
my $min_traces;       # what is the minimum number of traces that each species needs to have to continue processing it
my $species_list;     # optionally specify a file which contains species (quicker than looking up via separate script)
my $prog;             # specify path to a program that will produce list of (eukaryotic) species
my $ignore_processed; # check to see what files have previously been processed and ignore those even if gzip files are present
my $ls;               # just list what you will be fetching without actually fetching it
my $one_file;         # just grab files with the specified number
my $one_species;      # just grab files from the specified species

GetOptions ("max_files:i"      => \$max_files,
			"debug"            => \$debug,
			"timeout:i"        => \$timeout,
			"sleep:i"          => \$sleep,
			"max_attempts:i"   => \$max_attempts,
			"min_traces:i"     => \$min_traces,
			"species_list:s"   => \$species_list,
			"prog:s"           => \$prog,
			"ignore_processed" => \$ignore_processed,
			"ls"	           => \$ls,
			"one_file=s"       => \$one_file,
			"one_species=s"    => \$one_species);

# a quick couple of sanity checks on command line options
die "Use either -one_species <species name> or -species_list <file of species names>, but not both\n" if ($species_list && $one_species);
die "-one_file option must specify a 3 digit number (use leading zeroes if necessary)\n" if ($one_file && ($one_file !~ m/^[0-9]{3}$/));

###############################
#
# Set some default values
#
################################

$max_files = 20    if (!$max_files);
$timeout = 180     if (!$timeout);
$sleep = 10        if (!$sleep);
$max_attempts = 5  if (!$max_attempts);
$min_traces = 1000 if (!$min_traces);
$prog = glob("~/Work/bin/find_eukaryotes_in_trace_archive.pl -min_traces $min_traces") if (!$prog);
if (!$debug){$debug = 0}
else{$debug = 1}


#############################################################
#
# Check for previously processed files
#
#############################################################

# this will be the file that the downstream script (parse_ftp_tracedb_data.pl script) will read
my $processed_file_name = "trace_archive_processed_files.txt";

# if -ignore option is being used we need to first get a list of previously processed files
my %previously_processed;
if($ignore_processed){
	open(IN,"<$processed_file_name") or die "Can't find $processed_file_name file\n";
	while(<IN>){
		chomp;
		my ($file) = split(/\s+/,$_);
		$previously_processed{$file} = 1;
	}
	close(IN);
}

##############################################################
#
# Read in or generate a list of species to process
#
#############################################################

my @taxa = get_species_names();



########################
# BASIC FTP settings
########################
my $host = "ftp.ncbi.nlm.nih.gov";
my $user = "anonymous";
my $password = "krbradnam\@ucdavis.edu";
my $root = "pub/TraceDB";

# keep track of how many species did not have an exact file name match on FTP site
my $missing_counter = 0;
my $species_counter = 0;

my $ftp;
my $current_dir = getcwd;


##############################################################
#
# Main loop: ftp files for each species
#
#############################################################

SPECIES: foreach my $species (@taxa){
	
	$species_counter++;
	
	# make a directory for species if necessary
	my $path;
	$path = "$current_dir/$species";
	system("mkdir $path") unless -e $path;
	chdir $path;

	print STDERR "Processing files for $species\n";

	$ftp = Net::FTP->new($host, Debug => $debug, Timeout => $timeout)  or die "Cannot connect to $host: $@\n",$ftp->message;
	$ftp->login($user,$password) or die "Cannot login ", $ftp->message, "\n";
	$ftp->binary;

	my $dir = "$root/$species";


	# if we just want to list files, can stop here
	if($ls){
		my @files = $ftp->ls("$dir/fasta.$species.[0-9]*.gz");
		foreach my $file (@files){
			$file =~ s/.*fasta/fasta/;
			print "$file\n";
		}
		next SPECIES;
	}

	# need to find out how many files are in directory, grab all FASTA files and use array index of last file
	# also keep count of how many species we don't get an exact name match for
	my (@fasta);
	unless(@fasta = $ftp->ls("$dir/fasta.$species.[0-9]*.gz")){
		print STDERR "MISSING SPECIES: $species\n";
		$missing_counter++;
		next SPECIES;	
	}
		
	my @file_indices;
	foreach my $file (shuffle(@fasta)){
		my ($index) = $file =~ m/.*(\d{3})\.gz/;
		push(@file_indices,$index);
	}

	# if -one_file option is being used, can replace array with that one value
	@file_indices = ($one_file) if ($one_file);
	
	my $file_counter = 0;
	
	FILE: foreach my $index (@file_indices){
		$file_counter++;

		# break out of loop if we have exceeded max number of pages
		if ($file_counter > $max_files){
			print STDERR "Maximum number of files ($max_files) has been exceeded, skipping to next species\n";
			last FILE;
		}	
		
		# grab files in pairs, fasta + clip (they should pair up) 
		my ($fasta_return) = get_files($dir,$species,$index,"fasta",1);
		print STDERR "get_files (FASTA) failed for $species $index\n" if (!$fasta_return);
		
		my ($clip_return) = get_files($dir,$species,$index,"clip",1);
		print STDERR "get_files (CLIP) failed for $species $index\n" if (!$clip_return);
		
	}
	# tidy up
	$ftp->quit or die "Can't quit FTP",$ftp->message;

	print STDERR "\n";
}

print STDERR "\n$missing_counter species (out of $species_counter) could not be found on FTP site, might be due to slight variations in species names\n\n" if ($missing_counter);

exit;



##############################################################
#
#
#  T H E   S U B R O U T I N E S
#
#
#############################################################



sub get_files{
	my ($dir,$species,$index,$type,$attempt) = @_;

	# note that some files exist on the ftp site with a 'qualityless' part to their file name, e.g. clip.drosophila_melanogaster.qualityless.004.gz
	# these are rare and non-standard so we will just ignore them

	# format file name	
	my $file = "$type.$species.$index.gz";
	
	# now check to see whether this file has been processed before (if -ignore_processed option is in use)
	if($ignore_processed && $previously_processed{$file}){
		print STDERR "$file has been processed before, skipping to next file\n";
		return(1); 
	}
	
	# now need to check that files for species actually exist on FTP site
	if(defined $ftp->size("$dir/$file")){

		# get size of file on FTP site
		my $size = $ftp->size("$dir/$file");
		
		# is file in local directory AND same size?
		if(-e $file && (-s $file == $size)){
			print STDERR "$file exists locally - skipping\n";
			return(1);
		}
		# is file in local directory but different size?
		elsif(-e $file && (-s $file != $size)){
			print STDERR "FILE INCOMPLETE: refetching $file number $index, attempt number $attempt\n";	
		}
		# if we get here must be a new file to download
		else{
			print STDERR "fetching $file number $index, attempt number $attempt\n";	
		}

		# attempt to get file and use eval statement to catch any timeouts
		eval{$ftp->get("$dir/$file") or die "Can't grab $file\n",$ftp->message};	
		if ($@ =~ /Timeout/){
		   # catching code goes here
			print STDERR "$@\n";
			return(0) if ($attempt > $max_attempts);
			print STDERR "Sleeping for $sleep seconds, and then retrying\n";
			sleep($sleep);
			get_files($dir,$species,$index,$type,++$attempt);
		}		
	}
	# or give up
	else{
		print STDERR "MISSING FILE: $file not present on FTP site\n";
		# can't remember why I am returning 1 here...there was a reason
		return(1);
	}
	
	# if we get here, then we should have downloaded a file sucessfully
	return(1);
}



sub get_species_names{

	# either use a supplied list of species names, use the -one_species option, or run another program to generate them
	my @taxa;
	if($species_list or $one_species){
		if($species_list){
			print STDERR "Fetching list of species from $species_list\n\n";

			open (IN,"<$species_list") or die "Can't find file specified by -species_file: $species_list\n";

			while(my $species = <IN>){
				# skip blank-ish likes
				next unless $species =~ m/\w+/;
				# get species name in correct format (should already be lower case)
				chomp($species);	
				$species = lc($species);
				$species =~ s/ /_/g;
				push(@taxa,$species);
			}
			close(IN);			
		}
		elsif($one_species){
			@taxa = ($one_species);			
		}
	}
	else{
		print STDERR "Fetching list of eukaryotes by using $prog\n\n";
		@taxa = `$prog` or die "Can't run $prog\n";	
	}
	return(@taxa);
}



# signal event handler in case of interrupts (Ctrl+C)
sub INT_handler {
	
	# print final statistic of how many bases were clipped
	my $date = `date`; 
	chomp($date);
	
	print "\n\nSCRIPT INTERRUPTED at $date\n";
	print "\n$missing_counter species (out of $species_counter) could not be found on FTP site, might be due to slight variations in species names\n\n";    
	exit(0);
}