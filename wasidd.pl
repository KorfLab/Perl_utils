#!/usr/bin/perl
use strict; use warnings; use sigtrap;
use threads; use Thread::Queue;
#use File::Temp;
use Getopt::Std;
our ($opt_n, $opt_c, $opt_w, $opt_s, $opt_t, $opt_e, $opt_v);

##########################
### options processing ###
##########################

getopts('ncw:s:t:e:v');
my $WINDOW  = 4500;
my $STEP    = 1500;
my $THREADS = 1;
my $ENERGY  = "";
my $VERBOSE = 0;

die "WASIDD - DNA duplex destablilization predictor\n 
usage: wasidd.pl [options] <sequence> <working directory>
options:               [default]
  -n  nearest neighbor [copoly]
  -c  circular DNA     [linear]
  -w <int> window size [$WINDOW]
  -s <int> step size   [$STEP]
  -t <int> threads     [$THREADS]   (number of CPUs is optimal)
  -e <file> energies   [off] (tabular format: <position> <energy>)
  -v  verbose output   [off]
" unless @ARGV == 2;

my @ASIDD_OPT;
push @ASIDD_OPT, "-n" if $opt_n;
push @ASIDD_OPT, "-c" if $opt_c;
$THREADS = $opt_t if $opt_t;
$WINDOW  = $opt_w if $opt_w;
$STEP    = $opt_s if $opt_s;
$ENERGY  = $opt_e if $opt_e;
$VERBOSE = $opt_v if $opt_v;
die "window must be >= step" if $STEP > $WINDOW;

################################
### set up working directory ### 
################################

# globals defined in this section
my $FILE = $ARGV[0];
my $DIR  = $ARGV[1];
system("mkdir $DIR") unless -d $DIR;
my $DNA = "$DIR/DNA.raw";                  # file used for retrieving fragments
my $LENGTH;                                # total length of sequenced
my $ID;                                    # FASTA identifier

# sequence must be FASTA without ambiguities
open(DNA, ">$DNA") or die;
open(FASTA, $FILE) or die;
my $header = <FASTA>;
if ($header !~ /^>(\S+)/) {die "sequence must be in FASTA format\n"}
$ID = $1;
$LENGTH = 0;
while (<FASTA>) {
	if (/^>/) {die "only single FASTA sequences are allowed\n"}
	chomp;
	next unless /\S/;
	print DNA uc $_;
	$LENGTH += length($_);
}
close FASTA;
close DNA;
print STDERR "File: $FILE\n", "Sequence: $ID\n", "Length: $LENGTH\n",
	"Window: $WINDOW\n", "Step: $STEP\n", "Threads: $THREADS\n", if $VERBOSE;

###################
### energy file ### stored sparse, not many expected values
###################

my %Energy;
if ($ENERGY) {
	open(IN, $ENERGY) or die;
	while (<IN>) {
		my ($pos, $val) = split;
		$Energy{$pos - 1} = $val;
	}
	close IN;
}

#####################
### setup threads ###
#####################

my @idx = 0 .. ($LENGTH - $WINDOW + $STEP)/$STEP; # each idx is a window/job_id
my $total_jobs = @idx;
print STDERR "Jobs: $total_jobs\n" if $VERBOSE;
my $Q = new Thread::Queue; # thread-safe list
$Q->enqueue(@idx);
my @worker; # pool of workers
for (my $i = 0; $i < $THREADS; $i++) {$worker[$i] = threads->create(\&worker, $Q)}

#########################
### worker subroutine ###
#########################

sub worker {
	my ($q) = @_;
	
	my $tid = threads->tid;
	
	while ($q->pending) {
	
		# get a job
		my $idx = $q->dequeue;
		print STDERR $idx+1, "/$total_jobs (thread $tid): " if $VERBOSE;
		
		# skip if the job has already been done (allow restarts)
		my $size = -s "$DIR/$idx.asidd";
		$size = 0 if not defined $size;
		if ($size > 1000) {
			print STDERR "previously completed\n" if $VERBOSE;
			next;
		} elsif ($size > 0) {
			print STDERR "problematic, size=$size\n" if $VERBOSE;
			die "please check file $DIR/$idx.asidd and delete to restart\n";
		} else {
			print STDERR " processing\n" if $VERBOSE;
		}
		
		# get a window of DNA
		open(DNA, $DNA) or die;
		seek(DNA, $idx * $STEP, 0);
		my $seq;
		my $length = read(DNA, $seq, $WINDOW); # length of DNA fragment
		close DNA;
		die "crap, read error in thread $tid" if not defined $length;
		
		# create default analysis if containing runs of 10 or more Ns
		if ($seq =~ /N{10}/) {
			print STDERR "runs of Ns found, munging\n";
			open(ASIDD, ">$DIR/$idx.asidd") or die;
			for (my $i = 0; $i < $length; $i++) {
				print ASIDD "0\n";
			}
			close ASIDD;
			next;
		}
		
		# create temp seq file
		my $seqfile = "$DIR/seq.$tid.fasta";
		open(DNA, ">$seqfile") or die;
		print DNA ">window-$idx thread=$tid\n";
		for (my $i = 0; $i < length($seq); $i+=50) {
			print DNA substr($seq, $i, 50), "\n";
		}
		close DNA;
		
		# create temp energy file if there is data for it
		my @energy;
		for (my $i = 0; $i < $length; $i++) {
			my $position = $i + $idx * $STEP;
			if (exists $Energy{$position}) {
				push @energy, "$i\t$Energy{$position}";
			}
		}
		my $energyfile = "$DIR/energy.$tid";
		if (@energy) {
			open(OUT, ">$energyfile") or die;
			foreach my $e (@energy) {print OUT $e, "\n"}
			close OUT;
		}
		my $e = @energy ? "-e $energyfile" : "";
		
		# run asidd (always do this nicely)
		my $outfile = "$DIR/out.$tid";
		system("nice asidd @ASIDD_OPT $e $seqfile $length 1>$outfile 2>>$DIR/log");
		
		# parse output and save a slimmer version
		open(ASIDD, ">$DIR/$idx.asidd") or die;
		open(OUT, $outfile) or die;
		while (<OUT>) {
			if (/^position/) {last}
			elsif (/^Exceeded/) {
				# error condition met, create fake output and leave
				for (my $i = 0; $i < $length; $i++) {
					print ASIDD "0\n";
				}
				while (<OUT>) {} # just flushing the rest of the output
			}
		}
		while (<OUT>) {
			my ($position, $p, $g) = split;
			last if $position == $length;
			printf ASIDD "%.2f\n", $g;
		}
		close OUT;
		close ASIDD;
		
		# clean up
		unlink $seqfile    if -e $seqfile;
		unlink $outfile    if -e $outfile;
		unlink $energyfile if -e $energyfile;
	}
}

#######################
### process threads ###
#######################

# collect threads
for (my $i = 0; $i < $THREADS; $i++) {$worker[$i]->join}


# process results
print STDERR "processing results\n" if $VERBOSE;

open(DNA, $DNA) or die;
my $dna = <DNA>;
close DNA;

my $MID_POINT = $WINDOW/2;
my $LAST_FILE = int(($LENGTH-$WINDOW+$STEP)/$STEP);
my @asidd;

print "ASIDD output for $ID ($FILE)\n";
for (my $i = 0; $i < length($dna); $i++) {
	my $max = int($i / $STEP);
	my $min = int(($i - $WINDOW) / $STEP);
	$min = 0 if $min < 0;
	$max = $LAST_FILE if $max > $LAST_FILE;
	my $min_dist = 1e6;
	my $file;
	for (my $j = $min; $j <= $max; $j++) {
		my $midpoint = $j * $STEP + $MID_POINT;
		my $distance = abs($i - $midpoint);
		if ($distance < $min_dist) {
			$min_dist = $distance;
			$file = $j;
		}
	}
	
	# read and delete as necessary
	if (not defined $asidd[$file]) {
		print STDERR "." if $VERBOSE;
		$asidd[$file] = parse_asidd($file);
	}
	if ($min > 0 and defined $asidd[$min -1]) {
		undef $asidd[$min -1];
	}
	
	# create output at this position
	my $origin = $file * $STEP;
	my $offset = $i - $origin;
	my $g = $asidd[$file][$offset]; #"$file : $offset";
	print substr($dna, $i, 1), "\t", $g, "\n";
}
print STDERR "done\n" if $VERBOSE;


sub parse_asidd {
	my ($idx) = @_;
	my @data;
	open(ASIDD, "$DIR/$idx.asidd") or die "$DIR/$idx.asidd not found";
	while (<ASIDD>) {
		chomp;
		push @data, $_;
	}
	return \@data;
}
