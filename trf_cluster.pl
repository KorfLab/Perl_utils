#!/usr/bin/perl
#
# trf_cluster.pl
#
# A script to take a sample of tandem repeats and cluster them using BLAST
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;
use BPlite;
use IK;
use Cwd;
use DataBrowser;
use Getopt::Long;
use Keith;
use List::Util qw(sum shuffle);

########################
# Command line options
########################

my $local_peaks; # how many clusters to keep when analyzing BLAST results based on local alignments
my $global_peaks; # how many clusters to keep when analyzing BLAST results based on global alignments
my $mass_threshold; # what fraction of the maximum centromere mass should we use as a mimimum for retaining clusters
my $alignment_threshold; # minimum number of offset bases allowed to consider something a global alignment
my $blast_score; # -s parameter to use when running qstack
my $copies; # calculate z-axis of graph based on copy number rather than tandem repeat mass
my $data; # keep the data file that is produced, it's deleted by default
my $trf; # trf input file
my $sample; # how many (randomly chosen) Tandem repeats to use for BLAST cluster analysis
my $min_cluster_size; # how many sequences do you need to have in a cluster for it to count?
my $x_min;
my $x_max;
my $y_min;
my $y_max;
my $z_min;
my $z_max;
my $help; # get help

GetOptions ("local_peaks=i"    => \$local_peaks,
			"global_peaks=i"   => \$global_peaks,
			"mass_threshold=f" => \$mass_threshold,
			"alignment_threshold=i" => \$alignment_threshold,
			"blast_score=f" => \$blast_score,
			"min_cluster_size=i" => \$min_cluster_size,
			"copies" => \$copies,	
			"data" => \$data,		
			"trf=s" => \$trf,
			"sample=i" => \$sample,
			"x_min=i" => \$x_min,
			"x_max=i" => \$x_max,			
			"y_min=i" => \$y_min,
			"y_max=i" => \$y_max,
			"z_min=i" => \$z_min,
			"z_max=i" => \$z_max,
			"help" => \$help			
);

################
# set defaults
################
$local_peaks = 10         if (!$local_peaks);
$global_peaks = 10        if (!$global_peaks);
$alignment_threshold = 10 if (!$alignment_threshold);
$mass_threshold = 0.1     if (!$mass_threshold);
$blast_score = 50         if (!$blast_score);
$x_min = 1                if (!$x_min); 
$x_max = 600              if (!$x_max); 
$y_min = 20               if (!$y_min); 
$y_max = 80               if (!$y_max); 
$sample = 200000          if (!$sample);
$min_cluster_size = 50    if (!$min_cluster_size);

#############################
# check command-line options
##############################

my $usage = "
usage: trf-grapher.pl <options>
  -trf <trf output file>
  -local_peaks <int> [$local_peaks]
  -global_peaks <int> [$global_peaks]
  -alignment_threshold <int> [$alignment_threshold]
  -blast_score <float> [$blast_score]
  -copies : use copy number for the z-axis, rather than just tandem repeat mass
  -data : keep the data file produced for the graphs (deleted by default)
  -sample <int> : number of randomly chosen trace reads to use for BLAST analysis
  -help : this help
\n";

die $usage if ($help);
die $usage if (!$trf);

print STDERR "\n# $0 started at ", `date`, "\n";

######################################
#
# 1) Get species name from file name 
#
######################################

my $species;
if ($trf =~ /^([a-z]+_[a-z]+)/i) {$species = "$1"} 
else {$species = $trf}
my $species_name = $species;
$species_name =~ s/_/ /g;
$species_name =~ s/^(\w)/\u$1/;



#####################################################
#
# 2) setup directory and file names, hashes etc.
#
#####################################################

my $current_dir = getcwd;
my $path = "$current_dir/BLAST_analysis_$sample";
# make directory unless it already exists
system("mkdir $path") unless -e $path;

# will need two sample files
my $sample_trf = "$path/sample.trf";
my $sample_fa = "$path/sample.fa";

# will want to keep track of all the TIs that we sample from the FASTA file
my %sampled_tis; 


##############################################################################
#
# 3) count trace reads, create a sample fasta file, and calculate global GC% 
#
##############################################################################

# now make a selection of $sample sequences and write to new file

my ($reads, $nucleotides, $gc) = sample_fasta();

# if there are fewer tandems than requested, reset $sample
if($reads < $sample){
	print STDERR "\nWARNING: FASTA files only contains $reads sequences. This is less than requested size of $sample\n\n";
	$sample = $reads;
}



##############################################################################
#
# 4) Use combined TRF file to extract tandem repeats that are in sampled reads
#
##############################################################################

my $trf_count = parse_trfs();
print STDERR "Found $trf_count tandem repeats from $trf\n";


################################################################################
#
# 5) create X2 sequence file
#
################################################################################

# The X2 file is because the trf output will contain the predicted units of tandem repeats. If the tandem array
# in any one trace read does not have an 'edge', then it is impossible to know what the true boundary of each
# tandem repeat unit is. E.g. trf might say ATCG ATCG ATCG is a tandem array, but this might also just be 
# TCGA TCGA TCGA etc. The X2 file ensures that when we BLAST tandem repeats against each other, we will ensure
# we can cluster orthologous tandem repeats, even if their structure (as predicted by trf) looks very different

my $X2 = "$sample_trf.x2";
chdir $path;	

print STDERR "Creating X2 sequence file\n";
open(OUT, ">$X2") or die "Can't create $X2 file\n";

my %id_to_total_repeat_length; # total repeat length defined by a single tandem repeat
my %id_to_copy_number;         # number of copies of each tandem repeat
my %id_to_repeat_unit_length;  # length of repeat unit
my %id_to_gc;                  # GC content of tandem repeat unit
my %clean_data;                # for any data that survives cleaning, key is x-y coordinate pair info, value is 1 or 0
my $trf_counter = 0;           # keep track of how many repeats are processed
my ($max_x, $max_y) = (0, 0);

open(IN, $sample_trf) or die "Can't open $sample_trf file\n";	
my $fasta = new FAlite(\*IN);

print STDERR "Processing TRF file... ";

while (my $entry = $fasta->nextEntry) {
	$trf_counter++;

	my ($id, $copy_number, $duplicates, $unit_length, $repeat_fraction, $parent, $gc) = get_plot_data($entry);

	# add value of duplicates (if present) to trf count
	$trf_counter += $duplicates;
	
	# set max values of x & y
	($max_x = $unit_length) if ($unit_length > $max_x);
	($max_y = $gc) if ($gc > $max_y);
	
	$id_to_total_repeat_length{$id} = int($unit_length * $copy_number);
	$id_to_repeat_unit_length{$id}  = $unit_length;
	$id_to_copy_number{$id}         = $copy_number;
	$id_to_gc{$id}                  = $gc;
	
	# print to X2 file
	print OUT $entry->def, "\n", $entry->seq, "\n", $entry->seq, "\n";
}
close IN;
close(OUT);

print STDERR "done\n";



############################
#
# 6) create blast databases 
#
############################
print STDERR "Creating blast database\n";
if(!-e "$sample_trf.xnt"){
	system("xdformat -n -I $sample_trf") == 0 or die "Can't run xdformat to crate $sample_trf BLAST database\n";	
}
else{
	print STDERR "$sample_trf.xnt BLAST database arlready created, ";
}
print STDERR "done\n";



#################################################
#
# 7) all vs. all blast comparison @ 75% identity 
#
#################################################

print STDERR "Running qstack... ";
my $qstack = "tandem-vs-x2.qstack.gz";

if (!-e $qstack) {
	system("qstack.pl -dgs $blast_score -x10000 -i75 -w13 $sample_trf $X2 | gzip -c > $qstack") == 0 or die "Couldn't run qstack.pl\n";
} 
else {
	print STDERR " qstack already run, ";
}
print STDERR "done\n";



##############################
#
# 8) clustering and graphing 
#
##############################
# The result of this step is to create two cluster objects for eaah set of the top X global and local clusters
my $global_cluster = cluster_n_graph($qstack, 'global', $global_peaks);
my $local_cluster  = cluster_n_graph($qstack, 'local',  $local_peaks);




##############################
#
# 9) Write summary file 
#
##############################

open(OUT, ">${species}_summary.txt") or die;
print OUT "Summary for $species_name tandem repeat analysis\n";
printf OUT "%d reads, %d tandem repeats, %d bp, %.4f GC\n", $reads, $trf_count, $nucleotides, $gc;
my ($global_output, $global_top_length) = repeat_table('global', $global_cluster, $nucleotides, $global_peaks);
my ($local_output, $local_top_length) = repeat_table('local', $local_cluster, $nucleotides, $local_peaks);
print OUT "$global_output";
print OUT "$local_output";

# print out additional warning if top clusters from global and local modes are different sizes
print OUT "WARNING: Length of top cluster differs between global and local alignment modes\n" if ($global_top_length != $local_top_length);
print OUT "\nGlobal alignment cluster nodes\n", `cat global.cluster.seqs`;
print OUT "\nLocal alignment cluster nodes\n",`cat local.cluster.seqs`;
close OUT;

print STDERR "\n# $0 finished at ", `date`, "\n";




###############################################################################
#
#
#  T H E   S U B R O U T I N E S
#
#
###############################################################################



sub repeat_table {
	my ($mode, $cluster, $total, $limit) = @_;

	#  keep track of how many clusters will be skipped due to low mass or low size
	my $rejected = 0;

	# first want to find the max centromere mass value
	my $max_mass = 0;
	foreach my $cluster (@$cluster) {
		my $mass   = $cluster->{mass};
		if ($mass > $max_mass){
			$max_mass = $mass;
		}
	}
	
	my $threshold_percent = "%". ($mass_threshold * 100);
	
	# now loop through clusters again and print summary of each cluster as long as the length (i.e. centromere mass) of that
	# cluster is within a certain fraction of the maximum length. I.e. if there is only 1 significant cluster then we probably 
	# 	don't want to see details of 9 tiny little clusters
	my $output .= "\nRepeats ($mode mode) - using mass threshold $threshold_percent\n";
	
	# for global output we will include mean and stdev of pairwise identities
	if($mode eq "global"){
		$output .= join("\t", 'Cluster', 'n', 'id#', 'length', 'gc', 'mean', 'stdev', 'depth', 'mass', 'frac') . "\n";		
	} else {
		$output .= join("\t", 'Cluster', 'n', 'id#', 'length', 'gc', 'depth', 'mass', 'frac') . "\n";
	}
	
	my $count = 0;
	my $top_cluster_length = 0;
	
	foreach my $cluster (@$cluster) {
		my $size   = $cluster->{size};
		# no point going any further if cluster is too small
		if($size < $min_cluster_size){
			$rejected++;
			next;
		}
	
		
		# equally, no point continuing if mass is too low
		my $mass   = $cluster->{mass};
		if(($mass / $max_mass) <  $mass_threshold){
			$rejected++;
			next;
		}
	
		my $id     = $cluster->{id};
		my $depth  = $cluster->{depth};
		my $match  = $cluster->{match};
		my $gc     = $id_to_gc{$id};
		my $unit   = $id_to_repeat_unit_length{$id};
		my $frac = sprintf("%f", $mass/$total);
		
		if($mode eq "global"){
			my $mean   = $cluster->{mean};
			my $stdev  = $cluster->{stdev};
			$output .= join("\t", $count, $size, $id, $id_to_repeat_unit_length{$id}, $id_to_gc{$id}, $mean, $stdev, int($depth), $mass, $frac) . "\n";				
		}
		else{
			$output .= join("\t", $count, $size, $id, $id_to_repeat_unit_length{$id}, $id_to_gc{$id}, int($depth), $mass, $frac) . "\n";			
		}
		
		# want to grab the length of the top cluster
		$top_cluster_length = $id_to_repeat_unit_length{$id} if ($count == 0);
		$count++;
	}
	
	$output .= "\n";
	
	# check whether clusters 1 & 2 are very different in size and warn if not
	if(@$cluster > 2){
		my $size1 = ${$cluster}[0]{mass};
		my $size2 = ${$cluster}[1]{mass};
		my $threshold = 0.5;
		my $t_percent = $threshold * 100;
		$output .= "WARNING: Top clusters are similar in size. Cluster #1 has a tandem repeat mass within $t_percent% of cluster #0\n" if (($size2/$size1) > $threshold);		
	}
	if ($rejected){
		$output .= "NOTE: $rejected clusters are not reported due to having a tandem repeat mass less than $threshold_percent of cluster #0,\n";
		$output .= "or the cluster contained less than $min_cluster_size sequences\n"; 
		
	}
	$output .= "\n";
	
	return ($output,$top_cluster_length);
}
	
# Want to make clusters of similar tandem repeats
sub cluster_n_graph {
	my ($blast, $mode, $PEAKS) = @_;

	# keep track of the maximum centromere mass of any cluster
	my $max_mass = 0;
	
	# parse hits from blast output file
	print STDERR "Parsing blast report in $mode mode... ";
	my %hit;
	
	# will store details of percentage identity for all pairwise comparisons in a hash
	my %identity;
	
	open(IN, "gunzip -c $blast |") or die;
	while (<IN>) {
		# grab a subset of BLAST output:
		# query id, subject id, query start + end, subject start + end
		my @f = split;
		my ($q, $s, $identity, $qb, $qe, $sb, $se) = ($f[0],$f[1],$f[10],$f[17],$f[18],$f[20],$f[21]);
		my ($qid) = $q =~ /tandem-(\d+)/;
		my ($sid) = $s =~ /tandem-(\d+)/;
	
		# add percentage identity to hash, use query and subject IDs as combined key
		# each pair might have multiple matches, so just consider first one we come across as this will be the highest score
		unless (exists($identity{"${qid}_vs_${sid}"})){
			$identity{"$qid-$sid"} = $identity unless ($qid eq $sid); # no point adding self matches
		}
		
		# global vs local mode
		if ($mode eq 'global') {
			my $qlen = abs($qe - $qb) + 1;
			my $slen = abs($se - $sb) + 1;
			
			# calculate the difference between the repeat unit length and the length of the query and subject matches
			# if the entirety of a tandem repeat matches the query and subject, the difference will be zero
			my $qdif = abs($id_to_repeat_unit_length{$qid} - $qlen);
			my $sdif = abs($id_to_repeat_unit_length{$sid} - $slen);
			
			# only want to keep details of the BLAST match if it is a true global alignment
			# so reject if $qdif or $sdif is greater than some small value (default $alignment_threshold = 10)
			# but want to be a bit tighter if repeats are short, so effectively make threshold = 10% of 
			# of tandem repeat length for when the repeat length is <= 100 nt, otherwise use $alignment_threshold
			
			if($id_to_repeat_unit_length{$qid} <= 100){
				next if ($qdif > ($id_to_repeat_unit_length{$qid} * 0.1));
				
			} elsif ($qdif > $alignment_threshold){
				next;
			}

			if($id_to_repeat_unit_length{$sid} <= 100){
				next if ($sdif > ($id_to_repeat_unit_length{$sid} * 0.1));
				
			} elsif ($sdif > $alignment_threshold){
				next;
			}
#			next if ($qdif > $alignment_threshold or $sdif > $alignment_threshold);
		}
		# if we reach this point, we simply store the query and subject IDs of BLAST matches that we want
		# to keep for further processing, 
		$hit{$qid}{$sid} = 1;
	}
	close IN;
	print STDERR "done\n";
	
	print STDERR "Initial size of %hit hash is ", scalar(keys(%hit))," elements \n";
	
	# define clusters
	print STDERR "Clustering...\n\n";
	my @cluster;
	print "\n";
	my $loop_counter = 0;
	while (%hit) {
		$loop_counter++;
		print STDERR "Processing cluster $loop_counter in $mode mode\n";
		# find the 'depth' (total number of tandem repeats) and centromere mass for each head node
		my %depth;
		my %mass;
		foreach my $qid (keys %hit) {
			my $w = 0;
			my $tandem_repeat_mass = 0;
#			print "Q=$qid\n";
			foreach my $sid (keys %{$hit{$qid}}) {
				# for every match to the query sequence, we calculate the number of copies of tandem repeat
				# in each subject sequence and add those to a weighting factor $w
				# this means that one query could match just one subject but that subject occurred 100 times in one read
				# ...and this would be better than having 10 matches which only had 5 copies each in the original traces
				$w += $id_to_copy_number{$sid};
				$tandem_repeat_mass += $id_to_total_repeat_length{$sid}; 				
				#print "\tS=$sid ($id_to_copy_number{$sid}) W=$w\n";
			}

			# will default to using a weighting factor based on total copy number of any repeat across all BLAST matches
			$depth{$qid} = $w;				
			$mass{$qid} = $tandem_repeat_mass;				

			$max_mass = $tandem_repeat_mass if ($tandem_repeat_mass > $max_mass);
		}
					
		# we will assume that the set of blast hits with the highest 'depth' (stored in %depth)
		# will be the best place to start forming a cluster		

		# First find 'top hit' ($top). This will be the blast query ID with the highest weight 
		# There might be more than one with the same weight, in which case it will be randomly chosen???
		my ($top_depth) = sort {$depth{$b} <=> $depth{$a}} keys %hit;	
		my ($top_mass)  = sort {$mass{$b}  <=> $mass{$a}}  keys %hit;	
		
		# Now find all matching hits to top hit (including self match) and add to %match hash
		my %match;
		print STDERR "Top hit by mass is:  Q:$top_mass, mass = $mass{$top_mass} nt, total subject hits = ",scalar(keys(%{$hit{$top_mass}})),"\n";
		print STDERR "Top hit by depth is: Q:$top_depth, weight = $depth{$top_depth}, total subject hits = ",scalar(keys(%{$hit{$top_depth}})),"\n";

		# now just use mass as top hit
		my $top = $top_mass;

		foreach my $child (keys %{$hit{$top}}) {$match{$child} = 1}
#		print "1) Size of %hit hash is now ", scalar(keys(%hit))," elements \n";
	
		# Now if we assume that all of these matches to the top hit are related tandem repeats
		# then we can delete all occurences of these matches as other query sequences from the %hit hash
		foreach my $id1 (keys %match) {
#			print "\tDeleting Q:$id1 from %hit\n";
			delete $hit{$id1};
		} # parents
#		print "2) Size of %hit hash is now ", scalar(keys(%hit))," elements \n\n";
		
		# now delete the occurrence of any matches to the top hit (including self-hits) as matches to other query sequences in %hit
		# this might be quite rare but will happen, it won't change the number of primary keys in %hit though
		foreach my $id1 (keys %hit) {
			foreach my $id2 (keys %{$hit{$id1}}) {
				#print "\t\tDeleting the following match to Q:$id1: subject $id2\n" if (exists $match{$id2});
				delete $hit{$id1}{$id2} if exists $match{$id2}; # children
			}
		}

		# check for orphan entries in %hit, where all children have been removed
		foreach my $id1 (keys %hit) {
			if(scalar(keys%{$hit{$id1}}) == 0){
				#print "No children for Q:$id1, removing hash key\n";
				delete $hit{$id1};
			}
		}
		
		print STDERR "Size of %hit hash is now ", scalar(keys(%hit))," elements \n\n";
			
		# save cluster
		my @match = keys %match;
		push @cluster, {
			id => $top,
			depth => $depth{$top},
			match => \@match,
			mass => $mass{$top},
			size => scalar(@match),
		};		

	}
	print STDERR "done\n\n";

	# Extract sequences of head nodes
	# head nodes are the top query IDs of the N-best (default 10) clusters
	print STDERR "Saving head nodes... ";
	my $head = "$mode.cluster.ids";
	open(OUT, ">$head") or die;
	for (my $i = 0; $i < @cluster; $i++) {
#		last if $i == $PEAKS;
		last if (($cluster[$i]{mass} / $max_mass) < $mass_threshold);  
		print OUT "tandem-", $cluster[$i]{id}, "\n";
	}
	close OUT;
	print STDERR "done\n";
	
	my $seqs = "$mode.cluster.seqs";
	# grab sequences out of BLAST database
	system("xdget -n -f $sample_trf $head > $seqs") == 0 or die "Couldn't run xdget\n";

	# get sequences of each cluster and plot individual clusters
	print STDERR "Processing clusters...";
	for (my $i = 0; $i < @cluster; $i++) {
		#last if $i == $PEAKS;
		last if (($cluster[$i]{mass} / $max_mass) < $mass_threshold);  
		
		print STDERR " $i ";
	
		my $idfile     = "$mode.cluster$i.ids";
		my $seqfile    = "$mode.cluster$i.seqs";
		my $datafile   = "$mode.cluster$i.data";
	
		# make id file for cluster $i. This is a list of the matches to the top hit for that cluster
		open(LIST, ">$idfile") or die "Can't write to $idfile $!";
		
		# will also store each pairwise identity in an array
		my @identities = ();
		my $no_match = 0;
#		$cluster[$i]{size} = @{$cluster[$i]{match}};
#		print "Size = $cluster[$i]{size}\n";

		for (my $j = 0; $j < @{$cluster[$i]{match}}; $j++) {
			print LIST "tandem-${$cluster[$i]{match}}[$j]\n";

			# can also now calculate average pairwise identity for current cluster
			for (my $k = $j+1; $k < @{$cluster[$i]{match}}; $k++) {				
				my $key = ${$cluster[$i]{match}}[$j] . "-" . ${$cluster[$i]{match}}[$k];
				# not all posssible pairwise comparisons will be in the blast output,  
				# i.e. A matched B, and B matches C, but A does not match C even though they are all in cluster
				$no_match++ if (($mode eq "global") && !exists($identity{$key}));
				push(@identities,$identity{$key}) if (($mode eq "global") && exists($identity{$key}));
			}
		}
		close LIST;
 
		# calculate average pairwise identity and standard deviation (only if in global mode and you have at least 3 sequences in cluster)
		if(($mode eq "global") && ($cluster[$i]{size} >= $min_cluster_size)){
			my $pairwise_comparisons = @identities;	
			my $average_identity = sprintf("%.2f", (sum @identities) / $pairwise_comparisons);
#			print STDERR "\n$i) n=$cluster[$i]{size} There were $no_match pairwise comparisons that were not in BLAST output (compared to $pairwise_comparisons that were)\n";
#			print STDERR "Average \%identity = $average_identity\n";				
			
			$cluster[$i]{mean}  = $average_identity;
			$cluster[$i]{stdev} = sprintf("%.2f",sqrt(sum(map {($_ - $average_identity) ** 2} @identities) / ($pairwise_comparisons-1)));
			
		}
		
		
		
		# get sequences of matching IDs and write to file
		system("xdget -n -f $sample_trf $idfile > $seqfile") == 0 or die;
	
		# can now remove the ID file as we no longer need it (and IDs will be available in FASTA header of seqs file)
		unlink($idfile) or die "Can't remove $idfile\n";
		
		# make plot data for this cluster
		my @mat;
		open(IN, $seqfile) or die;
		my $fasta = new FAlite(\*IN);
		while (my $entry = $fasta->nextEntry) {
			my ($id, $copy_number, $duplicates, $unit_length, $repeat_fraction, $parent, $gc) = get_plot_data($entry);
			my ($x, $y, $z) = ($unit_length,$gc,$copy_number);			
			# weight $z to make it tandem repeat mass
			($z *= $x) unless ($copies);
			
			$mat[$x][$y] += $z;
		}
		close IN;
		
		# transform to match the aggregate data 
		# looks like we're just adding zero values on the z-axis rather than transforming??? [krb]
		for (my $x = 0; $x <= $max_x; $x++) {
			for (my $y = 0; $y <= $max_y; $y++) {

				if (not defined $mat[$x][$y]){
					$mat[$x][$y] = 0;			
				}
				else{
					# will change $z unless -copies is being used, effectively weight the copy number by the unit_length
					# and then calculate this mass as a fraction of the total amount of sequence processed
					# this will be a pseudo-approximation for what fraction of the genome is covered by tandem repeat
					unless ($copies){
						my $percent_tandem = ($mat[$x][$y]/$nucleotides) * 100;
						$mat[$x][$y] = $percent_tandem;
					}					
				}

			}
		}

		# write data file for possible graph for each cluster
		write_data_file($datafile, \@mat);	
	}

	# should now have lots of datafiles written to disk

	my $file = "aggregate.$mode";
#	plot("$file.script", "$file.pdf", @datafile);
	my $aggregate_data_file = "$file.data";
	# try making combined data file
	system("cat $mode.cluster[0-9]*.data > $aggregate_data_file") && die "Can't create combined data file\n";

	plot("$file.script", "$file.pdf", "$aggregate_data_file");
	unless ($data){
		foreach my $file (glob("*.data")){
			unlink("$file") or die "Can't remove $file\n";
		}
	}
	
	print STDERR "done\n\n";
	
	return \@cluster;
}

sub plot {
	my ($script, $output, $data_file) = @_;
	
	# form title for graph
	my $gbp = sprintf("%.2f",$nucleotides / 1000000000);
	my $title = "$species_name: $trf_counter repeats from $gbp Gbp.";
	 
	# process script file
	open(OUT, ">$script") or die;

	print OUT "library(scatterplot3d)\n";
	print OUT "pdf(\"$path/$output\")\n";
	print OUT "data <- read.table(\"$path/$data_file\")\n";
	print OUT "scatterplot3d(data,xlab=\"Repeat unit length (nt)\",ylab=\"GC%\",";
	if($copies){
		print OUT "zlab=\"Copy number\",";		
	}
	else{
		print OUT "zlab=\"%Tandem repeat\",";
	}

	print OUT "xlim=c($x_min,$x_max)," if ($x_min && $x_max);
	print OUT "ylim=c($y_min,$y_max)," if ($y_min && $y_max);
	print OUT "zlim=c($z_min,$z_max)," if ($z_min && $z_max);
	
	print OUT "highlight.3d=TRUE,col.axis=\"blue\",col.grid=\"lightblue\",main=\"$title\",type=\"h\",pch=20)\n";
	print OUT "dev.off()\n";
	close(OUT);
	
	# send output to /dev/null
	system("R CMD BATCH --vanilla $path/$script /dev/null") && die "Couldn't run \"R CMD BATCH --vanilla $path/$script /dev/null\"\n";
	
	# clean up files
	unlink("$script") or die "Can't remove $script\n";
	
}

###########################################################
#
# Grab a random selection of fasta sequences
#
###########################################################

sub sample_fasta{
	print STDERR "Processing all *processed_traces.nnn.fa files in current directory to sample $sample reads from\n";

	my $sample_nucleotides = 0;
	my $sample_reads = 0;
	my $sample_gc = 0;
	
	# first check to see whether a sample.fa already exists, if so we will use that and skip most of this step
	if (-e $sample_fa){
		print STDERR "Using existing $sample_fa file\n";
		open(IN, $sample_fa) or die "Can't open $sample_fa\n";
		my $fasta = new FAlite(\*IN);

		while (my $entry = $fasta->nextEntry) {
			my $seq = uc $entry->seq;
			my $header = $entry->def;
			my ($ti) = $header =~ m/ti\|(\d+) /;
			# keep the TI for later on
			$sampled_tis{$ti} = 1;
			
			# calculate some nucleotide stats for later on 
			$sample_reads++;
			$sample_nucleotides += length($seq);
			my $G = $seq =~ tr/G/G/;
			my $C = $seq =~ tr/C/C/;
			$sample_gc += $G + $C;
		}
		close IN;
		my $sample_genome_gc = $sample_gc/$sample_nucleotides;
		return($sample_reads,$sample_nucleotides, $sample_genome_gc);
	}
	

	# otherwise we have to loop through all fasta files


	# add all sequences to hashes (hope there is enough memory for this)
	my @fasta_files = glob("*processed_traces.[0-9][0-9][0-9].fa");
	
	# calculate how many sequences we need to extract from each file
	
	my $number_of_files = @fasta_files;
	my $seqs_per_file = int($sample / $number_of_files);

	print STDERR "Need to choose $seqs_per_file sequences from each trace read file\n";
	
	# two hashes to contain all of the sampled sequences
	my %ti_to_sampled_header;
	my %ti_to_sampled_seq;

	my $file_count = 0;
	my $selected = 0;
	
	foreach my $fasta_file (@fasta_files){
		$file_count++;
		print STDERR "\tProcessing $fasta_file\n";
		my %ti_to_header;
		my %ti_to_seq;

		open(IN, $fasta_file) or die "Can't open $fasta_file\n";
		my $fasta = new FAlite(\*IN);

		while (my $entry = $fasta->nextEntry) {
			my $header = $entry->def;
			my $seq = uc $entry->seq;
			my ($ti) = $header =~ m/ti\|(\d+) /;
			$ti_to_header{$ti} = $header;
			$ti_to_seq{$ti} = $seq;			
		}
		close IN;
		
		# shuffle TIs and choose random selection of $seqs_per_file sequence
		my @keys = shuffle(keys (%ti_to_header));
		my $number_of_tis = @keys;
		my $end = $seqs_per_file;

		# use another end point if we are at the last file in order to guarantee that we get all $sample sequences
		# we can have discrepancies as the number of seqs needed per file might be something like 150.4 and we can't get 0.4 
		# of a sequence. So when we process the last file, we can take a few more sequences if necessary
		$end = ($sample - $selected) if ($file_count == $number_of_files);

		# If the file does not contain enough sequences, we will just all sequences in that file
		$end = $number_of_tis if ($seqs_per_file > $number_of_tis);
		
		# now select the TIs of our chosen number
		my @selected = @keys[0..$end-1];
		
		# @selected now contains all of the TIs that we need from the current file, loop through these and add to hash
		foreach my $ti (@selected){
			$selected++;
			$ti_to_sampled_header{$ti} = $ti_to_header{$ti};
			$ti_to_sampled_seq{$ti} = $ti_to_seq{$ti};
		}
	}
	
	# at this point %ti_to_sampled_header, and %ti_to_sampled_seq should contain all we need to make our sample file
	# so loop through these and print them out and also add details to $nucleotides, $reads, and $gc
	
	print STDERR "Creating new FASTA file of sampled processed trace reads\n";
	open(OUT,">$sample_fa") or die "Can't create $sample_fa\n";
	
	foreach my $ti (keys (%ti_to_sampled_seq)){
		my $header = $ti_to_sampled_header{$ti};
		my ($ti) = $header =~ m/ti\|(\d+) /;
		# keep the TI for later on
		$sampled_tis{$ti} = 1;

		my $seq = uc($ti_to_sampled_seq{$ti});
		my $tidied_seq = Keith::tidy_seq($seq);
		print OUT "$header\n$tidied_seq\n";

		# calculate some nucleotide stats for later on 
		$sample_reads++;
		$sample_nucleotides += length($seq);
		my $G = $seq =~ tr/G/G/;
		my $C = $seq =~ tr/C/C/;
		$sample_gc += $G + $C;
	}
	close OUT;
	my $sample_genome_gc = $sample_gc/$sample_nucleotides;
	return($sample_reads, $sample_nucleotides, $sample_genome_gc);
}

###########################################################
#
# Grab a selection of TRFs from the combined file
# This is to speed up the BLAST step
#
###########################################################

sub parse_trfs{

	print STDERR "Creating sample TRF file\n";

	my $trf_counter = 0;
	
	# write TRF sequences to new file and note the TI number in another hash
	open(OUT,">$sample_trf") or die "Can't create $sample_trf\n";
	open(IN, $trf) or die "Can't open $trf file\n";	
	my $fasta = new FAlite(\*IN);
	while (my $entry = $fasta->nextEntry) {
		my $header = $entry->def;
		my ($ti) = $header =~ m/ti\|(\d+) /;

		# was this TI in our sample of trace reads? If so, print out tandem repeat ot sample TRF file
		if($sampled_tis{$ti}){
			$trf_counter++;
			my $seq = $entry->seq;
			my $tidied_seq = Keith::tidy_seq($seq);
			print OUT "$header\n$tidied_seq\n";
		}
	}
	close IN;
	close OUT;
	return($trf_counter);
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
		# just use a zero value for duplicates field when using a non-slimmed TRF file (which won't have duplicate info)
		my ($id, $copy_number, $unit_length, $repeat_fraction, $parent) = $entry->def =~ />tandem-(\d+) N=(\S+) L=(\d+) F=(\d+)% P=(\S+)/;
		return ($id, $copy_number, 0, $unit_length, $repeat_fraction, $parent, $gc);            
    }
}


sub write_data_file {
	my ($file, $m) = @_;
	open(OUT, ">$file") or die;
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
