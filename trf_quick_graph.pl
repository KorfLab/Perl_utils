#!/usr/bin/perl
#
# trf_quick_graph.pl
#
# A script to make an x,y,z plot from a set of Tandem Repeats
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use strict;
use warnings;

use FAlite;
use BPlite;
use IK;
use Cwd;
use DataBrowser;
use Getopt::Long;

########################
# Command line options
########################

my $trf; # trf input file
my $species; # specify species name (optional)
my $clean; # remove final counts of repeat based on %percentage of the most frequent repeat
my $copies; # calculate z-axis of graph based on total number of repeat copies, rather than tandem repeat mass
my $peakdump; # dump a FASTA file of all sequences under the highest peak
my $data; # keep the data file that is produced, it's deleted by default
my $x_min;
my $x_max;
my $y_min;
my $y_max;
my $help; # get help

GetOptions ("trf=s" => \$trf,
			"species=s" => \$species,
			"clean" => \$clean,
			"copies" => \$copies,
			"x_min=i" => \$x_min,
			"x_max=i" => \$x_max,			
			"y_min=i" => \$y_min,
			"y_max=i" => \$y_max,
			"peakdump" => \$peakdump,
			"help" => \$help,
			"data" => \$data			
);


#############################
# check command-line options
##############################

my $usage = "
usage: trf-grapher.pl <options>
  -trf <trf output file>
  -species <name of species>
  -clean : use a 1% theshold to clean up data
  -copies : calculate z-axis as total tandem repeat copy number (default is to use tandem repeat mass)
  -data : keep data file used for graph
  -peakdump : dump out sequences from top peak
  -x_min, -x_max, -y_min, -y_max <set axes ranges>
  -help - this help
\n";

die $usage if ($help);
die $usage if (!$trf);
$species = "speciesX" if (!$species);
$x_min = 1               if (!$x_min); 
$x_max = 600             if (!$x_max); 
$y_min = 20              if (!$y_min); 
$y_max = 80              if (!$y_max);

#############################
# setup directory and files #
#############################

my $path = getcwd;


###################################
# 1) read repeats from trf output
###################################

my @matrix;
my ($max_x, $max_y, $max_z) = (0, 0, 0);

my $trf_counter = 0;

open(IN, $trf) or die "Can't open $trf file\n";
my $fasta = new FAlite(\*IN);

while (my $entry = $fasta->nextEntry) {
	$trf_counter++;
	my ($id, $copy_number, $duplicates, $unit_length, $repeat_fraction, $parent, $gc) = get_plot_data($entry);
	
	# also add duplicates to $trf counter
	$trf_counter += $duplicates;
	
	# add info to a 3-dimensional matrix
	my ($x, $y, $z) = ($unit_length, $gc, $copy_number);
	
	# unless -copies is being used will change $z by effectively weighting the copy number by the unit_length
	unless ($copies){
		$z *= $x;
	}
	
	$matrix[$x][$y] += $z;
	
	# set max values of x & y
	($max_x = $x) if ($x > $max_x);
	($max_y = $y) if ($y > $max_y);
}
close IN;
 


########################
# 2) Fill @matrix array
########################
my ($x_maxima, $y_maxima) = (0,0);
# All axes will have some undefined values, e.g. there may not be any repeats at 45 nt length
# or repeats at 50 nt length will not possess all possible GC% contents
# so loop through all *possible* x and y values and then add zero counts to corresponding z value
for (my $x = 0; $x <= $max_x; $x++) {
	for (my $y = 0; $y <= $max_y; $y++) {
		$matrix[$x][$y] = 0 if not defined $matrix[$x][$y];
		if ($matrix[$x][$y] > $max_z){
			$max_z = $matrix[$x][$y];
			($x_maxima,$y_maxima) = ($x,$y);
		}
	}
}

print STDERR "Peak is at $x_maxima nt, $y_maxima GC%\n";

###############################################
# 3) create aggregate data file and make plot
###############################################

# if using -clean option clean things up by removing low copy repeats from @matrix
if($clean){
	# take 1% of the maximum z-axis as threshold value
	my $clean_threshold  = int($max_z * 0.01);						
	
	print "Using -clean option to remove any data with a Z-axis value less than $clean_threshold\n";
	for (my $x = 0; $x <= $max_x; $x++) {
		for (my $y = 0; $y <= $max_y; $y++) {
			$matrix[$x][$y] = 0 if not defined $matrix[$x][$y];
			$matrix[$x][$y] = 0 if ($matrix[$x][$y] < $clean_threshold);
		}
	}
}

# this step uses R to do the plotting. First writes a script which is sent to R
write_data_file("${trf}.data", \@matrix);
plot("${trf}.script", "${trf}.pdf", "${trf}.data");


# do we want to print out the sequences that contribute to the highest peak?
if($peakdump){
	open(IN, $trf) or die "Can't open $trf file\n";
	my $fasta = new FAlite(\*IN);
	open(OUT, ">$trf.peak") or die "Can't create $trf.peak file\n";
	while (my $entry = $fasta->nextEntry) {
		my ($id, $copy_number, $duplicates, $unit_length, $repeat_fraction, $parent, $gc) = get_plot_data($entry);
		# now we just want to skip to the next sequence unless the x & y coords correspond to maxima values
		next unless ($unit_length == $x_maxima && $gc == $y_maxima);
		print OUT $entry->def, "\n", $entry->seq, "\n";					
	}
	close IN;
	close OUT;
}

exit;

###############################################################################
# subroutines
###############################################################################

sub plot {
	my ($script, $output, $data_file) = @_;

	# calculate sequence size in Gbp
	my $title = "$species: $trf_counter repeats. Maxima at $x_maxima nt, $y_maxima% GC";
	$title .= ", CLEAN" if ($clean);

	open(OUT, ">$script") or die "Couldn't open $script file\n";
	print OUT "library(scatterplot3d)\n";
	print OUT "pdf(\"$path/$output\")\n";
	print OUT "data <- read.table(\"$path/$data_file\")\n";
	print OUT "scatterplot3d(data,xlab=\"Repeat unit length (nt)\",ylab=\"GC%\",";
	if($copies){
		print OUT "zlab=\"Tandem Repeat Copy number\",";		
	}
	else{
		print OUT "zlab=\"Tandem Repeat Mass\",";
	}

	print OUT "xlim=c($x_min,$x_max)," if ($x_min && $x_max);
	print OUT "ylim=c($y_min,$y_max)," if ($y_min && $y_max);
	
	print OUT "highlight.3d=TRUE,col.axis=\"blue\",col.grid=\"lightblue\",main=\"$title\",type=\"h\",pch=20)\n";
	print OUT "dev.off()\n";
	close(OUT);
	
	# send output to /dev/null
	system("R CMD BATCH --vanilla $path/$script /dev/null") && die "Couldn't run \"R CMD BATCH --vanilla $path/$script /dev/null\"\n";

	# clean up files
	unlink("$script") or die "Can't remove $script\n";
	(unlink("$data_file") or die "Can't remove $data_file\n") unless ($data);

}

# get basic tandem repeat information from trf_wrapper.pl output file
sub get_plot_data {
	my ($entry) = @_;
	my $seq = uc $entry->seq;
	my $gc = gc($seq);

	# need to handle the data differently if it comes from the the 'slim' version of the TRF output file 
	if($trf =~ m/high.slim.trf/){
		my ($id, $copy_number, $duplicates, $unit_length, $repeat_fraction, $parent) = $entry->def =~ />tandem-(\d+) N=(\S+) D=(\d+) L=(\d+) F=(\d+)% P=(\S+)/;
		return ($id, $copy_number, $duplicates, $unit_length, $repeat_fraction, $parent, $gc);
	}
	else{
		# just use a zero value for duplicates field
		my ($id, $copy_number, $unit_length, $repeat_fraction, $parent) = $entry->def =~ />tandem-(\d+) N=(\S+) L=(\d+) F=(\d+)% P=(\S+)/;
		return ($id, $copy_number, 0, $unit_length, $repeat_fraction, $parent, $gc);		
	}
}


sub write_data_file {
	my ($file, $m) = @_;
	open(OUT, ">$file") or die "Can't write to $file\n";
	for (my $x = 0; $x < @$m; $x++) {
		for (my $y = 0; $y < @{$m->[$x]}; $y++) {
#			print OUT "$x $y $m->[$x][$y]\n";
			# print out x,y,z data, but don't print out values if Z-axis (copy number) is zero
			print OUT "$x $y $m->[$x][$y]\n" unless ($m->[$x][$y] == 0);
		}
	}
	close OUT;
}
sub gc {
	my ($seq) = @_;
	my $A = $seq =~ tr/A/A/;
	my $C = $seq =~ tr/C/C/;
	my $G = $seq =~ tr/G/G/;
	my $T = $seq =~ tr/T/T/;
	return int(100 * ($G+$C) / ($A+$C+$G+$T));
}
