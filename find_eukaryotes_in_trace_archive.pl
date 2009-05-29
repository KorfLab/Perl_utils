#!/usr/bin/perl
#
# find_eukaryotes_in_trace_archive.pl
#
# A script to produce a list of all species in the NCBI trace archive, and then cross reference them 
# using info from two dumped files from the NCBI Taxonomy database, and only report on those which are eukaryotes
#
# Last updated by: $Author$
# Last updated on: $Date$


use strict; 
use warnings; 
use Getopt::Long;

my $verbose;    # turn on extra output - e.g. errors for traces that were rejected
my $program;    # can specify path to query_tracedb.pl program
my $min_traces; # minimum number of traces that need to be present in a species to consider it worthy of mentioning

GetOptions ("verbose"      => \$verbose,
            "program:s"    => \$program,
			"min_traces:i" => \$min_traces);


# set defaults if not specified on command line
$min_traces = 1000 if (!$min_traces);
$program = "query_tracedb.pl" if (!$program);


# 1st hash will contain species names from names.dmp (key = species name, value = tax ID)
# 2nd hash will be the reverse look up (key = tax ID, value = species name)
# 3rd hash will only contain those names which we know to be eukaryotes
# 4th hash will track synonyms, as query_tracedb & the NCBI FTP site don't always use the strict taxonomic name for a species
# this makes an assumption that synonyms are unique names
my %species_name_to_id;
my %id_to_species_name;
my %eukaryotes;
my %synonym_to_id;

# expecting to find the following two files in the current directory (from NCBI taxonomy dump file)
my ($fnames, $fnodes) = ("names.dmp", "nodes.dmp");

# parse names file
open(IN, $fnames) or die "$fnames not found";
while (<IN>) {
	# keep everything lower case
	$_ = lc($_);
	
	# only want to process lines which deal with scientific names of a species or it's synonyms
	next unless (/scientific name/ || /synonym/);
	my ($taxid, $name) = split(/\|/, $_);
	$name =~ s/^\s+//g;
	$name =~ s/\s+$//g;
	$taxid =~ s/\s//g;

	# we either have a proper name for a species (in which case put details in two hashes)
	if (m/scientific name/){
		$species_name_to_id{$name} = $taxid;
		$id_to_species_name{$taxid} = $name;
	}
	# or we have a synonym
	else{
		$synonym_to_id{$name} = $taxid;
	}
}
close IN;


# parse nodes file
my (%Parent, %Child);
open(IN, $fnodes) or die;
while (<IN>) {
	$_ = lc($_);
	my ($taxid, $parent) = split(/\|/, $_);
	$taxid =~ s/\s//g;
	$parent =~ s/\s//g;
	if ($taxid == $parent) {
		next if ($taxid == 1);
		die "error child == parent $taxid";
	}
	$Parent{$taxid} = $parent;
	$Child{$parent}{$taxid} = 1;
}
close IN;

# create relationship tree
my $root = 1;  # eukaryota = tax ID 2759
my $tree = createTree(\%Parent, \%Child, $root, {});

# now descend into tree from the root of eukaryotes (2759) to find names of eukaryotes
descendTree($tree->{1}{131567}{2759});


# use the stat function of query_tracedb to get a list of all species represented in trace archive
# this will give species names + trace count for that species
my @tracedb_output = `$program stat`;

my $species_counter = 0;
my $eukaryote_counter = 0;
my $min_traces_counter = 0;

SPECIES: foreach my $species (@tracedb_output){
	$species_counter++;
	chomp($species);
	$species = lc($species);

	# skip the first and last lines of output
	next if ($species =~ m/^organism/);
	next if ($species =~ m/total:/);

	# now need to separate sequence count from species name
	$species =~ s/\s+([0-9,]+)$//;
	my $count = $1;
	die "$species\n" if (!$count);
	$count =~ s/,//g;

	if($eukaryotes{$species}){
		# keep track if there are too few traces, and skip to next species
		if ($count < $min_traces){
			$min_traces_counter++; 
			next SPECIES;
		}
		$eukaryote_counter++;
		
		if($verbose){
			print "$species ($count traces) is a EUKARYOTE\n";			
		}
		else{
			print "$species\n";
		}
	
	}
	# if species name is not in main hash, need to look up in synonyms hash
	elsif(exists($synonym_to_id{$species})){
			# keep track if there are too few traces, and skip to next species
			if ($count < $min_traces){
				$min_traces_counter++; 
				next SPECIES;
			}
			$eukaryote_counter++;
			if($verbose){
				print "$species is a synonym of \'$id_to_species_name{$synonym_to_id{$species}}\' ($count traces), which is a EUKARYOTE\n";				
			}
			else{
				print "$species\n";
			}
	}
	else{
		print "$species - ($count traces) is NOT a eukaryote!\n" if ($verbose);
	}
}

if($verbose){
	print "\n\n";
	print "Found $eukaryote_counter eukaryotes out of $species_counter species present in NCBI trace archive\n";
	print "$min_traces_counter eukaryotes were ignored for having fewer than $min_traces traces\n";
}
exit;


sub descendTree {
	my ($tree) = @_;
	foreach my $child (keys %$tree) {
#		($eukaryotes{lc($id_to_species_name{$child})} = 1) if scalar keys %{$tree->{$child}} == 0;
		($eukaryotes{lc($id_to_species_name{$child})} = 1);
		descendTree($tree->{$child});
	}
}


# recursive funciton for taxonomy tree
sub createTree {
	my ($P, $C, $parent, $tree) = @_;
	foreach my $child (keys %{$C->{$parent}}) {
		$tree->{$parent}{$child} = {};
		createTree($P, $C, $child, $tree->{$parent});
	}
	return $tree;
}

