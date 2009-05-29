#!/usr/bin/perl
#
# parse_ftp_tracedb_data.pl
#
# A script to 
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Keith;
use FAlite;
use Getopt::Long;
use Net::FTP;

# turn on autoflushing off buffer
$| = 1;

# Need this or else FTP won't work!
$ENV{FTP_PASSIVE} = 1;

$SIG{'INT'} = 'INT_handler'; # for handling signal interrupts


my $dir;              # directory to look for *.gz files
my $min_bases;        # minimum number of bases that you need in a sequence after clipping to keep it
my $max_n;	          # what is the maximum percentage of N's allowed in the clipped sequence
my $max_output_size;  # maximum size of output files (in gigabytes) before they are rotated
my $ignore_processed; # check to see what files have previously been processed and ignore those even if gzip files are present
my $verbose;          # turn on extra output - e.g. errors for traces that were rejected
my $move_files;       # ftp files to commando
my $clean_files;      # remove existing gzip files after transfer
my $stop;         # stop script when you reach species starting with specified letters

GetOptions ("dir:s"             => \$dir,
			"min_bases:i"       => \$min_bases,
			"max_n:f"           => \$max_n,
			"verbose"           => \$verbose,
			"max_output_size:f" => \$max_output_size,
			"ignore_processed"  => \$ignore_processed,
			"move_files"        => \$move_files,
			"stop:s"            => \$stop,                     
			"clean_files"       => \$clean_files);


# set defaults if not specified on command line
$min_bases = 30    if (!$min_bases);
$max_n = 5         if (!$max_n);
die "-stop option must specify a lower case letter\n" if($stop && ($stop !~ m/[a-z]/));

# set an initial limit of 0.5 Gb, and will create new output files if they grow bigger than that
$max_output_size = 0.5 if (!$max_output_size);
# convert to bytes
$max_output_size *= 1073741824;

$dir = "" if (!$dir);

# add a trailing slash if none was specified to path (this is a bit of a kludge)
($dir .= "/") if ($dir && $dir !~ m/\/$/);


# need to keep track of:
# 1) total traces parsed
# 2) how many bases are clipped
# 3) how many traces are rejected for being too short after clipping low quality bases, 
# 4) how many traces are rejected for containing too many Ns
# 5) how many other errors there were (e.g. when clip left coordinate is greater than total length of sequence)
my $total_traces = 0;
my $total_bases = 0;
my $total_clipped_bases = 0;
my $total_rejected_high_n = 0;
my $total_rejected_too_short = 0;
my $total_errors = 0;

my $date;


################################################################################
#
#                          M A I N   P R O G R A M
#
################################################################################

# open a file which will contain details of any files that have been processed
# (so they can be ignored in future)
my $processed_file_name = "trace_archive_processed_files.txt";


if(-e "$processed_file_name"){
	open(PROCESSED,">>$processed_file_name") or die "Can't append to $processed_file_name\n";
}
else{
	open(PROCESSED,">$processed_file_name") or die "Can't create $processed_file_name\n";
}

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

# Get list of clip files in current directory (or in $dir if specified)
my @clip_files = glob("${dir}clip.*.gz");

FILE: foreach my $clip_file (@clip_files) {
	
	# make copy of just file name without path
	my $clip_file_name = $clip_file;
	$clip_file_name =~ s/^$dir//;
	
	# extract species name and file number
	$clip_file_name =~ m/clip.(.*?)\..*?([0-9]+).gz/;
	my ($species,$file_number) = ($1,$2);
	
	# stop if -stop is being used
    if($stop && ($species =~ m/^$stop/)){
    	print "\nLetter $stop has been reached...stopping program\n";
        exit(0);
     }
	
	# make versions for fasta file name
	my ($fasta_file, $fasta_file_name) = ($clip_file, $clip_file_name);
	$fasta_file =~ s/clip/fasta/;
	$fasta_file_name =~ s/clip/fasta/;

	
	# no point proceeding if we can't find the accompanying FASTA file in the same directory
	if (! -e $fasta_file){
		print "ERROR: Can't find $fasta_file\n"; 
		next FILE;
	}

	# now check to see whether this file has been processed before (if -ignore_processed option is in use)
	if($ignore_processed && $previously_processed{$clip_file_name}){
		print "$clip_file_name has been processed before, skipping to next file\n" if ($verbose);
		next FILE;	
	}

  	print "Unzipping $clip_file\n" if ($verbose);	

	# now process clip file to store information
	# use two hashes which track clip left/right values for each TI
	my %ti_to_clip_left;
	my %ti_to_clip_right;

	open(IN, "gunzip -c $clip_file |") || die "Can't open $clip_file: $? $!\n";	
	while(my $line = <IN>){
		next if ($line =~ m/^TI/);
		my ($ti,$left,$right) = split (/\s+/,$line);
		$ti_to_clip_left{$ti} = $left;
		$ti_to_clip_right{$ti} = $right;
	}
	close(IN);

	print "Unzipping $fasta_file\n" if ($verbose);	
  
	# now need to parse the associated data in FASTA file
	# send command through a pipe to gunzip and then output will be sent to FAlite module	
	# want to store fasta headers, fasta sequences, and seq lengths in hashes all tied to TI number of each sequence            
	my %ti_to_seq;
    my %ti_to_header;
    my %ti_to_seqlength;

	# reset file-level counters
	my $file_total_bases = 0;
	my $file_total_traces = 0;
	my $file_clipped_bases = 0;
	my $file_rejected_high_n = 0;
	my $file_rejected_too_short = 0;
	my $file_errors = 0;

    open(FASTA,"/usr/bin/gunzip -c $fasta_file | ") or die "Can't open pipe on gunzip $fasta_file: $? $!\n";
    my $FA = new FAlite (\*FASTA);
    while (my $entry = $FA->nextEntry) {
            my $header = $entry->def;
            $header =~ m/ti\|(\d+) /;
            my $ti = $1;
            die "No ti in $header\n" if (!$ti);
            $ti_to_header{$ti} = $header;
            $ti_to_seq{$ti} = $entry->seq;
            $ti_to_seqlength{$ti} = length($entry->seq);

			#print "TI: $ti CLIP_RIGHT: $ti_to_clip_right{$ti}\n";
			
            # update stats
			$total_traces++;
			$file_total_traces++;
            $total_bases += length($entry->seq);
			$file_total_bases += length($entry->seq);
    }       
    close(FASTA);


	# time to open an output file for processed data, but we first need to see whether there are output files 
	# that were already copied to commando and/or output files in the current directory which were not 
	# copied to commando. This will help decide what numerical suffix the file should get
	my $last_commando_file = check_commando($species);
	my $last_local_file = check_local($species,$dir);
	
	# choose highest value between two files
	my $last_processed_output_file = $last_commando_file+1; 
	($last_processed_output_file = $last_local_file) if ($last_local_file > $last_commando_file); 

	# increment file counter (it will be set to 1 if there are no files on commando)
	my $output_file_counter = $last_processed_output_file;
	my $output_file = "${species}_processed_traces.${output_file_counter}.fa";
	
	# if output file doesn't already exist open a new one...
	if(! -e $output_file){
		open(OUT,">$output_file") or die "Can't create $output_file\n";		
	}
	# so we already have an existing output file for this species, but do we want to use it?
	else{
		# if output file exists but is too big we need to make a new file with an incremented suffix
		if(-s $output_file > $max_output_size){
			print "$output_file has become too large (>$max_output_size gigabytes), creating new output file\n" if ($verbose);
			
			# do we want to move this older file to commando?
			if ($move_files){
				# first zip file to save space and speed transfer
				system("/usr/bin/gzip $output_file") && die "Could not gzip $output_file\n";
				ftp_files("${output_file}.gz");
				
				# if we are cleaning then we can also get rid of the processed gzip file
				if($clean_files){
					unlink("${output_file}.gz")  or die "Can't remove ${output_file}.gz\n";
				}
				
			}
		
			$output_file_counter++;
			$output_file = "${species}_processed_traces.${output_file_counter}.fa";
			open(OUT,">$output_file") or die "Can't write to $output_file\n";		
		}
		# if here then we have an existing output file which is not too big, so we can just append 
		else{
			open(OUT,">>$output_file") or die "Can't append to $output_file\n";			
		}
	}

	
	print "Processing $clip_file_name & $fasta_file_name\n";	

	# now loop through all the TIs and clip sequence if necessary
	TRACE: foreach my $ti (sort {$a <=> $b} keys %ti_to_seq){
		my $length = $ti_to_seqlength{$ti};

		# set clip left value to zero and set clip right-values to the length of sequence as we want to see if there are values that are lower than this		
		my $left=0;
		my $right=$length;
		
		my $clip_left = $ti_to_clip_left{$ti};
		my $clip_right = $ti_to_clip_right{$ti};

		# now assign $left the most extreme values discovered
		($left = $clip_left)     if ($clip_left > $left);

		# some clip_right values are greater than the length of the sequence, in which case we can change to
		# set them to the sequence length, i.e. no right clipping
		($clip_right = $length) if ($clip_right > $length);

		# for the right-fields, need to remember that a zero value just means no clip information is present
		($right = $clip_right)   if (($clip_right < $right)   && ($clip_right   != 0));

#		print "TI=$ti Length = $length Left = $left, right = $right\n";


		# have some basic sanity checks to catch errors		
		if(($clip_left > $length) || ($left > $right)){
			$file_errors++;
			$total_errors++;
			print "ERROR: Inconsistent information - TI:*$ti* LENGTH: $length CLIP_LEFT:$clip_left CLIP_RIGHT:$clip_right\n" if ($verbose);
		
			# no point going any further
			next TRACE;
		}
		
		# now clip sequence if necessary
		my $seq = $ti_to_seq{$ti};
		if($left > 0 || $right < $length){
			
			# want to keep track total number of clipped bases and species specific clipped bases
			my $clipped_bases = $left + ($length - $right - 1);

			$total_clipped_bases += $clipped_bases;
			$file_clipped_bases += $clipped_bases;

			# how many bases will be left after clipping?
			my $remaining_bases = $length - $clipped_bases;
			
			# Add check in case remaining sequence is below some useful limit?		
			if($remaining_bases < $min_bases){
				$total_rejected_too_short++;
				$file_rejected_too_short++;
				print "ERROR: $ti has $remaining_bases bases after clipping\n" if ($verbose);
				next(TRACE);
			}
			
			# use substr to do the actual clipping and modify corresponding FASTA header with clipping info
			$seq = substr($seq,$left-1,$remaining_bases);
			$ti_to_header{$ti} .= " CLIPPED: $clipped_bases nt";
		}

		# die if we don't have any sequence for some reason
		die "No seq\n$ti\tlength=$length\tleft=$left\tright=$right\n" if (!$seq);

		# now check whether remaining sequence contains more than $max_n Ns in which case we should ignore it
		my $n = ($seq =~ tr/N/N/);
		my $percent_n =($n / length($seq)*100);
		if($percent_n > $max_n){
			$total_rejected_high_n++;
			$file_rejected_high_n++;			
			print "ERROR: $ti contains more than ${max_n}% Ns in its sequence\n" if ($verbose);
			next(TRACE);
		}
		
		# If we have survived to this point, then we have a valid sequence which we can tidy and then print to output
		my $tidied_seq = tidy_seq($seq);
		print OUT "$ti_to_header{$ti}\n$tidied_seq\n";
	}
	close(OUT);
	
	# print summary statistics for just this file
	my $percent_clipped = sprintf("%.1f",($file_clipped_bases/$file_total_bases)*100);
	print "Processed $file_total_traces traces containing $file_total_bases nt of which $file_clipped_bases nt (${percent_clipped}%) had to be clipped\n";
	print "$file_rejected_too_short traces were rejected for being too short (<$min_bases) after vectory/quality clipping\n" if ($file_rejected_too_short > 0);
	print "$file_rejected_high_n traces were rejected for containing too many unknown bases (>$max_n%)\n" if ($file_rejected_high_n > 0);
	print "$file_errors traces were rejected for containing errors (inconsistant information)\n" if ($file_errors > 0);

	# update processed file information, include settings used to process files
	print PROCESSED "$clip_file_name min_bases=$min_bases max_n=$max_n\n";
	print PROCESSED "$fasta_file_name min_bases=$min_bases max_n=$max_n\n";
	
	
	# check to see if there are more files to come this species, if not then potentially FTP the output file to commando
	# first have to form what the next file name will be (for clip files)
	my $next_number = $file_number;
	$next_number =~ s/^0+//;
	$next_number++;
	$next_number = sprintf("%03d", $next_number);
	my $next_file_name = "clip.$species.$next_number.gz";

	# transfer file now unless there is another file for this species in the same directory
	unless(-e $next_file_name){
		system("/usr/bin/gzip $output_file") && die "Could not gzip $output_file\n";
		if ($move_files){
			ftp_files("${output_file}.gz");
			
			# if we are cleaning then we can also get rid of the processed gzip file
			if($clean_files){
				unlink("${output_file}.gz")  or die "Can't remove ${output_file}.gz\n";
			}						
		}
	}
	
	print "\n";
	
	# do we want to remove the original (unprocessed) zipped files?
	if ($clean_files){
		unlink($clip_file)  or die "Can't remove $clip_file\n";
		unlink($fasta_file) or die "Can't remove $clip_file\n";		
	}
	
}

close(PROCESSED);

# print summary statistics for all files examined so far, across all species
my $percent_clipped = sprintf("%.1f",($total_clipped_bases/$total_bases)*100);
print "\n\n======================================================\n\n";
print "TOTAL: Processed $total_traces traces containing $total_bases nt of which $total_clipped_bases nt (${percent_clipped}%) had to be clipped\n";
print "TOTAL: $total_rejected_too_short traces were rejected for being too short (<$min_bases) after vectory/quality clipping\n";
print "TOTAL: $total_rejected_high_n traces were rejected for containing too many unknown bases (>$max_n%)\n";
print "TOTAL: $total_errors traces were rejected for containing errors (inconsistant information)\n\n";
print "======================================================\n\n";

exit(0);




# signal event handler in case of interrupts (Ctrl+C)
sub INT_handler {
	
	# print final statistic of how many bases were clipped
	$date = `date`; 
	chomp($date);
	
	my $percent_clipped = sprintf("%.1f",($total_clipped_bases/$total_bases)*100);
	print "\n\n======================================================\n\n";
	print "SCRIPT INTERRUPTED at $date\n";
	print "TOTAL: Processed $total_traces traces containing $total_bases nt of which $total_clipped_bases nt (${percent_clipped}%) had to be clipped\n";
	print "TOTAL: $total_rejected_too_short traces were rejected for being too short (<$min_bases) after vectory/quality clipping\n";
	print "TOTAL: $total_rejected_high_n traces were rejected for containing too many unknown bases (>$max_n%)\n";
	print "TOTAL: $total_errors traces were rejected for containing errors (inconsistant information)\n\n";
	print "======================================================\n\n";
    exit(0);
}

sub ftp_files{
	
	my @files = @_;

	########################
	# BASIC FTP settings
	########################
	my $host = "commando.genomecenter.ucdavis.edu";
	my $user = "tracedb";
	my $password = "Korf2009";
	my $dir = "tracedb_lite";

 	my $timeout = 180;
	my $ftp;

	$ftp = Net::FTP->new($host, Timeout => $timeout)  or die "Cannot connect to $host: $@\n",$ftp->message;
	$ftp->login($user,$password) or die "Cannot login ", $ftp->message, "\n";
	$ftp->binary;
	$ftp->cwd($dir) or die "Can't change directory to $dir",$ftp->message;

	# now transfer all files
	foreach my $file (@files){
		print "FTPing $file to $host\n";
		$ftp->put($file) or die "Can't transfer $file to $host",$ftp->message;	
	}

	$ftp->quit or die "Can't quit FTP",$ftp->message;

}

sub check_commando{
	
	my $species = shift;

	########################
	# BASIC FTP settings
	########################
	my $host = "commando.genomecenter.ucdavis.edu";
	my $user = "tracedb";
	my $password = "Korf2009";
	my $dir = "tracedb_lite";

 	my $timeout = 180;
	my $ftp;

	$ftp = Net::FTP->new($host, Timeout => $timeout)  or die "Cannot connect to $host: $@\n",$ftp->message;
	$ftp->login($user,$password) or die "Cannot login ", $ftp->message, "\n";
	$ftp->binary;
	$ftp->cwd($dir) or die "Can't change directory to $dir",$ftp->message;

	my @fasta;
	# if there no files already copied to commando, we return zero
 	unless(@fasta = $ftp->ls("${species}_processed_traces.[0-9]*.fa.*")){
		$ftp->quit or die "Can't quit FTP",$ftp->message;
		return(0);               
 	}               

	# else we grab the numerical suffix of the last file there and return it's suffix instead
    my $last_file = $fasta[-1];
    $last_file =~ m/${species}_processed_traces.([0-9]*).fa.*/;
	my $suffix = $1;
	$ftp->quit or die "Can't quit FTP",$ftp->message;
    return($suffix);

}

sub check_local{
	
	my $species = shift;
	my $dir = shift;

	my @fasta;
	# if there no files already copied to commando, we return zero
 	unless(@fasta = glob("${dir}${species}_processed_traces.[0-9]*.fa")){
		return(0);               
 	}               

	# else we grab the numerical suffix of the last file there and return it's suffix instead
    my $last_file = $fasta[-1];
    $last_file =~ m/${species}_processed_traces.([0-9]*).fa/;
	my $suffix = $1;
    return($suffix);

}

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