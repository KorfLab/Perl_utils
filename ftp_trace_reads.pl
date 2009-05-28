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

# turn on autoflush
$| = 1;

# for handling signal interrupts
$SIG{'INT'} = 'INT_handler'; 

# Need this or else things won't work!
$ENV{FTP_PASSIVE} = 1;

my $max_files;    # maximum number of files to download for any species
my $debug;        # whether to turn on debugging in FTP module
my $timeout;      # set timeout value for Net::FTP module
my $sleep;        # how long to sleep for before retrying ftp download
my $max_attempts; # how many attempts to download one file before giving up
my $min_traces;   # what is the minimum number of traces that each species needs to have to continue processing it
my $species_list; # optionally specify a file which contains species (quicker than looking up via separate script)
my $prog;         # specify path to a program that will produce list of (eukaryotic) species
my $verbose;      # turn on extra output (usefulf for debugging)
my $ignore_processed; # check to see what files have previously been processed and ignore those even if gzip files are present
my $stop;         # stop script when you reach species starting with specified letters

GetOptions ("max_files:i"    => \$max_files,
			"debug"          => \$debug,
			"timeout:i"      => \$timeout,
			"sleep:i"        => \$sleep,
			"max_attempts:i" => \$max_attempts,
			"min_traces:i"   => \$min_traces,
			"species_list:s" => \$species_list,
			"stop:s"         => \$stop,			
			"prog:s"         => \$prog,
			"ignore_processed"  => \$ignore_processed,			
			"verbose"        => \$verbose);

# set defaults if not specified on command line
$max_files = 2     if (!$max_files);
$timeout = 180     if (!$timeout);
$sleep = 10        if (!$sleep);
$max_attempts = 5  if (!$max_attempts);
$min_traces = 1000 if (!$min_traces);
die "-stop option must specify a lower case letter\n" if($stop && ($stop !~ m/[a-z]/));
$prog = glob("~/Work/bin/find_eukaryotes_in_trace_archive.pl -min_traces $min_traces") if (!$prog);

if (!$debug){
	$debug = 0;
}
else{
	$debug = 1;
}


#########################################
# GET LIST OF SPECIES TO IGNORE
#########################################

# this will be the file that parse_ftp_tracedb_data.pl will produce
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

#########################################
# GET LIST OF EUKARYOTES IN TRACE ARCHIVE
#########################################

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

SPECIES: foreach my $species (@taxa){
	
	# stop if -stop is being used
	if($stop && ($species =~ m/^$stop/)){
		print "Letter $stop has been reached...stopping program\n";
		exit(0);
	}
	$species_counter++;
	
	print "\nFetching files for $species\n" if ($verbose);

	$ftp = Net::FTP->new($host, Debug => $debug, Timeout => $timeout)  or die "Cannot connect to $host: $@\n",$ftp->message;
	$ftp->login($user,$password) or die "Cannot login ", $ftp->message, "\n";
	$ftp->binary;

	my $dir = "$root/$species";

	# need to find out how many files are in directory, grab all FASTA files and use array index of last file
	# also keep count of how many species we don't get an exact name match for
	my @fasta;
	unless(@fasta = $ftp->ls("$dir/fasta.$species.[0-9]*.gz")){
		print "MISSING SPECIES: $species\n";
		$missing_counter++;
		next SPECIES;	
	}		
		
	# need to work out what is the first and last files in the directory 
	# for some species the first file is not numbered 1 (it's either absent or there are qualityless files instead of
	# the regular files)
	
	my $first_file = $fasta[0];	
	$first_file =~ m/fasta.$species.(\d+).gz/;
	my $starting_file = $1;
	$starting_file =~ s/^0+//;
	
	
	my $last_file = $fasta[-1];
	$last_file =~ m/fasta.$species.(\d+).gz/;
	my $number_of_files = $1;
	
	$number_of_files =~ s/^0+//;
	
	my $file_counter = 0;
	
	FILE: for (my $i=$starting_file;$i<=$number_of_files;$i++){
		
		$file_counter++;
						
		# break out of loop if we have exceeded max number of pages
		if ($file_counter > $max_files){
			print "Maximum number of files ($max_files) has been exceeded, skipping to next species\n" if ($verbose);
			last FILE;
		}	
		
		# grab files in pairs, fasta + clip (they should pair up) 
		my $fasta_return = get_files($dir,$species,$i,$number_of_files,"fasta",1);
		print "get_files (FASTA) failed for $species $i\n" if (!$fasta_return);
		
		my $clip_return = get_files($dir,$species,$i,$number_of_files,"clip",1);
		print "get_files (CLIP) failed for $species $i\n" if (!$clip_return);
		
	}
	# tidy up
	$ftp->quit or die "Can't quit FTP",$ftp->message;
}



print "\n$missing_counter species (out of $species_counter) could not be found on FTP site, might be due to slight variations in species names\n\n";

exit;



sub get_files{
	my ($dir,$species,$i,$number_of_files,$type,$attempt) = @_;

	# note that some files exist on the ftp site with a 'qualityless' part to their file name, e.g. clip.drosophila_melanogaster.qualityless.004.gz
	# these are rare and non-standard so we will just ignore them

	# format file name	
	my $number = sprintf("%03d", $i);
	my $file = "$type.$species.$number.gz";

	# now check to see whether this file has been processed before (if -ignore_processed option is in use)
	if($ignore_processed && $previously_processed{$file}){
		print "$file has been processed before, skipping to next file\n" if ($verbose);
		return(1);
	}

	
#	print "Looking for $file\n";
	
	# now need to check that files for species actually exist on FTP site
	if(defined $ftp->size("$dir/$file")){

		# get size of file on FTP site
		my $size = $ftp->size("$dir/$file");
		
		# is file in local directory AND same size?
		if(-e $file && (-s $file == $size)){
			print "$file exists locally - skipping\n\n" if ($verbose);
			return(1);
		}
		# is file in local directory but different size?
		elsif(-e $file && (-s $file != $size)){
			print "FILE INCOMPLETE: refetching $file number $i of $number_of_files, attempt number $attempt\n";	
		}
		# if we get here must be a new file to download
		else{
			print "fetching $file number $i of $number_of_files, attempt number $attempt\n";	
		}

		# use eval statement to catch any timeouts
		eval{$ftp->get("$dir/$file") or die "Can't grab $file\n",$ftp->message};	
		if ($@ =~ /Timeout/){
		   # catching code goes here
			print "$@\n";
			return(0) if ($attempt > $max_attempts);
			print "Sleeping for $sleep seconds, and then retrying\n";
			sleep($sleep);
			get_files($dir,$species,$i,$number_of_files,$type,++$attempt);
		}		
	}
	
	# or give up
	else{
		print "MISSING FILE: $file  not present on FTP site\n";
		# can still return 1 as this i
		return(1);
	}
	

	# if we get here, then we should have downloaded a file sucessfully
	return(1);
}



sub get_species_names{
	print "\nFetching list of eukaryotes by using $prog\n\n";

	my @taxa;
	if($species_list){
		open (IN,"<$species_list") or die "Can't find file specified by -species_file: $species_list\n";
		while(my $species = <IN>){
			
		# get species name in correct format (should already be lower case)
		chomp($species);	
		$species = lc($species);
		$species =~ s/ /_/g;
		push(@taxa,$species);
		}
		close(IN);
	}
	else{
		@taxa = `$prog`or die "Can't run $prog\n";	
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