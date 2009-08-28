#!/usr/bin/perl
#
# This script finds any overlapping genes in a worm genome, then test those genes to find
# if genes overlap, or if CDS overlap, or if transcript overlaps.  Also finds
# exons associated with each of the genes. Then compares the exons to determing if
# genes are located wholly within an intron or if genes are possibly interdigitating.
#
# Written by Yi Zhang, July/August 2009
#
use strict; 
use warnings; 
use Getopt::Long; 
use worm_GFFlite;

# user must specify a directory containing valid GFF files, other options can be used
die "usage: $0 \noptions: \n --overwrite => create and overwrite subsetfiles \n --lines => prints out lines of output \n --stats => print out statistics \n --threshold => threshould value\n --web => give weblinks\n --output => output data to file\n --rna => output all data <GFF directory>\n\n" unless @ARGV > 0;

my %gene_to_start = ();	
my %gene_to_stop = (); 
my %gene_to_strand = (); 	
my %gene_to_CDSs = ();	# key is gene id and value is a list of CDS ids belonging to that gene

my %CDS_to_transcript = ();	# key is CDS id and value is a list of transcript ids belonging to that CDS
my %CDS_to_gene = ();
my %CDS_to_start = (); 	
my %CDS_to_stop = (); 	
my %CDS_to_exons = ();	# key is CDS id, value is list of exon coordinates belonging to the CDS
my %CDS_to_introns = ();	# key is CDS id, value is list ofintron coordinates belonging to the CDS
my %CDS_to_strand = ();

my %exon_to_start = ();	
my %exon_to_stop = ();
my %intron_to_start = ();
my %intron_to_stop = ();

my %transcript_to_start = ();
my %transcript_to_stop = ();
my %transcript_to_strand = ();

### counters for all the different types of overlap ###
my $counter_gene_overlap = 0;
my $gene_same = 0;
my $gene_diff = 0;
my $counter_rna_overlap = 0;
my $counter_rna_rna_overlap = 0;
my $rna_rna_same_strand = 0;
my $rna_rna_diff_strand = 0;
my $counter_rna_proteingene_overlap = 0;
my $rna_proteingene_same_strand = 0;
my $rna_proteingene_diff_strand = 0;
my $counter_protein_coding_overlap = 0;
my $counter_gene = 0;
my $counter_transcript_CDS_overlap = 0;
my $transcript_CDS_same = 0;
my $transcript_CDS_diff = 0;
my $counter_transcript_overlap = 0;
my $transcript_same = 0;
my $transcript_diff = 0;
my $counter_gene_in_intron = 0;
my $intron_same = 0;
my $intron_diff = 0;
my $counter_interdigitates = 0;
my $interdigitates_same = 0;
my $interdigitates_diff = 0;

########## 1. COMMAND LINE OPTIONS ##########
my $threshold;	# establish a threshold, default threshold distance should be about 50000
my $web;	# option to output weblinks for the genes
my $output;	# allows user to output the result to textfile.
my $rna;	# prints all outputs (including gene comparisons) if option is turned on
my $overwrite;	# creates new subsetfile that will overwrite the previous subsetfile
my $stats; # prints out the stats
my $lines; # prints out lines of output


GetOptions (
	"threshold=i" => \$threshold,
	"web!" => \$web,
	"output:s" => \$output,
	"rna!" => \$rna,
	"overwrite!" => \$overwrite,
	"stats!" => \$stats,
	"lines!" => \$lines,
); 

die "Please use either --lines or --stats on the command line\n" if (!$stats && !$lines); # no stats or line option, then no output
die "--web option is used with --lines option\n" if ($web && !$lines); # weblinks modifies the output lines

$threshold = 50000 if (!$threshold); # set threshold default at 50,000
if ($web) {my $base_url = "http://www.wormbase.org/db/gene/gene?name=" }; # set up url if option's turned on
# user has option to output to file or to screen
# if user outputs to file, gives user the option to specify file name to output to
# otherwise it will output to default file named output.txt
if (defined ($output)) {
	if (!$output) { 
		my $output = "output.txt";
		open (OUT, ">$output") or die "Can not output to file\n"; 
	}	 	
	if ($output) {
		open (OUT, ">$output") or die "Can not output to file\n"; 
	}		
}

######## 2.a. ACCESSING DIRECTORY ##########
my $dir = $ARGV[0];
my @filtered_files; # example: subfile_CHROMOSOME_I.gff
opendir (DIR, "$dir") or die "can not open directory $dir\n";
if ($overwrite) {
	while ($_ = readdir(DIR)) {
		if ($_ =~ /\S+\.gff/ and $_ !~ /subfile_(\S+\.gff)/) { 
			print "overwriting file $_...";
			my $overwritten_file = &GFFlite::gff_reader_alt ($dir , $_);
			push (@filtered_files, $overwritten_file);
			print "finished overwriting\n";
		}
	}
	close DIR; print "\n";
}
my @unfiltered_files = (); # example: CHROMOSOME_I.gff
if (!$overwrite){
	while ($_ = readdir(DIR)) {
		if ($_ =~ /.*?\.gff/) {
			if ($_ =~ /subfile_(\S+\.gff)/) { 
				push (@filtered_files, "$dir"."$_"); #list of all sub_files: dir/subfile_file.gff
			}
			else {push (@unfiltered_files, "$_");} #list of all unprocessed gff files: dir/file.gff
		}
	}
	close DIR; 
	foreach my $var (@unfiltered_files) {
		my $var_filtered = "$dir"."subfile_"."$var";
		if (-e $var_filtered){next;}
		else {
			print "filtering file $var...\n";
			my ($dir_output) = GFFlite::gff_reader_alt ($dir , $var);
			push (@filtered_files, $dir_output);
			print "\tfile $var_filtered created\n";
		}
	}
}

# prints column headers on screen or to file
print "TYPE\tLENGTH_OF_OVERLAP\tGENE1\tTRANSCRIPT1\tGENE2\ttTRANSCRIPT2\tSTRANDS\tCOMMENT\n" if (!defined ($output) && !$web && $lines);	
print "TYPE\tLENGTH_OF_OVERLAP\tSTRANDS\tCOMMENT\tURL\n" if (!defined ($output) && $web && $lines);
print OUT "TYPE\tLENGTH_OF_OVERLAP\tGENE1\ttTRANSCRIPT1\ttTRANSCRIPT2\tCDS2\tSTRANDS\tCOMMENT\n" if (defined ($output) && !$web && $lines);
print OUT "TYPE\tLENGTH_OF_OVERLAP\tSTRANDS\tCOMMENT\tURL\n" if (defined ($output) && $web && $lines);

########## 2.c. PROCESS FILES ##########
# read the filtered files in directory
foreach my $gff_subfile (@filtered_files) {	
	print "Processing file: $gff_subfile\n" if (defined ($output) && $lines);
	# put into subroutine to generate hashes
	&file_reader_alt ($gff_subfile);
	print "\tFile processed and closed, comparing genes begins...\n" if (defined ($output) && $lines) ;
	$gff_subfile =~ /.*?subfile_(\S+)\.gff/; my $chromosome = $1;
	########## 3. TESTING FOR OVERLAP ##########
	# create key to look up what's returned
	# note: since it has been sorted by start value, options 2 and 3 are no longer possible
	my %overlap_description = (
		0 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tno overlap between geneA and geneB",
		1 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\t3'end of geneA overlaps 5'end of geneB",
		2 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\t5'end of geneA overlaps 3' end of geneB",
		3 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tgeneB contains geneA",
		4 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tgeneA contains geneB",
		5 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tgeneA and geneB are identical",	
		6 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tgeneA contained in intron of geneB", #CDS/Intron
		7 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tgeneB contained in intron of geneA", #Intron/CDS
		8 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\texons overlap", 
		9 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tgeneA and geneB possibly interdigitates", #Gene/Gene
		10 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tno exon overlap, geneA and geneB interdigitates", 
	);
	
	my @genelist; # empty list to store genes after they're sorted
	
	# sorts the hash by increasing numeric start value, and stores the key to genelist
	foreach my $key (sort {$gene_to_start{$a} <=> $gene_to_start{$b}} keys %gene_to_start) {push @genelist, $key;}
	
	OUTTER: while (@genelist) {
		
		my $geneA = shift (@genelist); # $geneA represents the 1st gene
	 
		INNER: foreach my $j (@genelist) { # $j represents the 2nd gene
		
			# make sure that the two genes to be tested are w/in the threshold
			next OUTTER if (abs($gene_to_start{$geneA}-$gene_to_start{$j})>$threshold);
			
			# passing 2 sets of gene start/stop coordinates to subroutine to check for overlap
			my ($a,$b) = &overlap_check($gene_to_start{$geneA}, $gene_to_stop{$geneA}, $gene_to_start{$j}, $gene_to_stop{$j});
			next INNER if ($a == 0);	# if there's no overlap, move to next $j gene
			
		# counters======================
			$counter_gene_overlap++;
			$gene_same++ if ($gene_to_strand{$geneA} eq $gene_to_strand{$j});
			$gene_diff++ if ($gene_to_strand{$geneA} ne $gene_to_strand{$j});
		#===============================
	
			# if there's any overlap but the genes doesn't have CDSs 
			# if rna option is turned on, prints out the type of overlap b4 moving to next gene
			if (!exists $gene_to_CDSs{$geneA} or !exists $gene_to_CDSs{$j}) {
				
			# counters======================
				$counter_rna_overlap++; 
				if (!exists $gene_to_CDSs{$geneA} and !exists $gene_to_CDSs{$j}) {
					$counter_rna_rna_overlap++;
					if ($gene_to_strand{$geneA} ne $gene_to_strand{$j}) {$rna_rna_same_strand++;}
					if ($gene_to_strand{$geneA} eq $gene_to_strand{$j}) {$rna_rna_diff_strand++;}
				}
				if (exists $gene_to_CDSs{$geneA} or exists $gene_to_CDSs{$j}) {
					$counter_rna_proteingene_overlap++;
					if ($gene_to_strand{$geneA} ne $gene_to_strand{$j}) {$rna_proteingene_same_strand++;}
					if ($gene_to_strand{$geneA} eq $gene_to_strand{$j}) {$rna_proteingene_diff_strand++;}
				}				
			#===============================
				
				# shows rna overlaps if option is turned on
				if ($rna) {	
					my $type;
					if (!exists $gene_to_CDSs{$geneA} && exists $gene_to_CDSs{$j}) {$type = "RNA/Gene"}
					if (exists $gene_to_CDSs{$geneA} && !exists $gene_to_CDSs{$j}) {$type = "Gene/RNA"}
					if (exists $gene_to_CDSs{$geneA} && exists $gene_to_CDSs{$j}) {$type = "RNA/RNA"}
					my $description1_return = &description1 ($overlap_description{$a}, $geneA, $j);
					print OUT "$type\t$b\t$description1_return\n" if (defined ($output));
					print "$type\t$b\t$description1_return\n" if (!defined ($output));
				}
				next;
			}
		
			$counter_protein_coding_overlap++ if (exists $gene_to_CDSs{$geneA} && exists $gene_to_CDSs{$j});
			
			# place the two genes in subroutine to check for CDS overlap, transcript overlap, etc...
			my $CDS_check_return = &CDS_check($geneA, $j);
			if ($CDS_check_return ne 'null' && $lines) {
				print OUT "$CDS_check_return\n" if (defined ($output));
				print "$CDS_check_return\n" if (!defined ($output));
			}
		}
	}
	print "\t\tComparing genes in file complete\n\n" if (defined $output && $lines);
}

# statistics are always printed to output whether or not --stats option is on.
if (defined ($output)) {
	print OUT "Total number of gene: $counter_gene\n";
	print OUT "Total number of gene overlaps: $counter_gene_overlap\n";
	print OUT "\tsame strand: $gene_same\n";
	print OUT "\tdifferent strand: $gene_diff\n";
	print OUT "Total number of rna gene overlap: $counter_rna_overlap\n";
	print OUT "\trna overlapping rna: $counter_rna_rna_overlap\n\t\tsame strand: $rna_rna_same_strand\n\t\tdifferent strand: $rna_rna_diff_strand\n";
	print OUT "\trna and protein coding gene overlaps: $counter_rna_proteingene_overlap\n\t\tsame strand: $rna_proteingene_same_strand\n\t\tdifferent strand: $rna_proteingene_diff_strand\n";
	print OUT "Total number of protein coding overlap: $counter_protein_coding_overlap\n";
	print OUT "\tTotal number of gene in intron: $counter_gene_in_intron\n\t\tsame strand: $intron_same\n\t\tdifferent strand: $intron_diff\n";
	print OUT "\tTotal number UTR, UTR overlap: $counter_transcript_overlap\n";
	print OUT "\t\tsame strand: $transcript_same\n\t\tdifferent strand: $transcript_diff\n";
	print OUT "\tTotal number UTR, CDS overlap: $counter_transcript_CDS_overlap\n";
	print OUT "\t\tsame strand: $transcript_CDS_same\n\t\tdifferent strand: $transcript_CDS_diff\n";
	print OUT "\tTotal possible interdigitating genes: $counter_interdigitates\n\t\tsame strand: $interdigitates_same\n\t\tdifferent strand: $interdigitates_diff\n\n";
}

if (defined ($output)) {
	if (!$output) { print "Output to output.txt complete\n\n";}	 	
	if ($output) {print "Output to $output complete\n\n";}		
	close (OUT);
}

# prints statistics out to screen if option is on
if ($stats) {
	print "Total number of gene: $counter_gene\n";
	print "Total number of gene overlaps: $counter_gene_overlap\n";
	print "\tsame strand: $gene_same\n";
	print "\tdifferent strand: $gene_diff\n";
	print "Total number of rna gene overlap: $counter_rna_overlap\n";
	print "\trna overlapping rna: $counter_rna_rna_overlap\n\t\tsame strand: $rna_rna_same_strand\n\t\tdifferent strand: $rna_rna_diff_strand\n";
	print "\trna and protein coding gene overlaps: $counter_rna_proteingene_overlap\n\t\tsame strand: $rna_proteingene_same_strand\n\t\tdifferent strand: $rna_proteingene_diff_strand\n";
	print "Total number of protein coding overlap: $counter_protein_coding_overlap\n";
	print "\tTotal number of gene in intron: $counter_gene_in_intron\n\t\tsame strand: $intron_same\n\t\tdifferent strand: $intron_diff\n";
	print "\tTotal number UTR, UTR overlap: $counter_transcript_overlap\n";
	print "\t\tsame strand: $transcript_same\n\t\tdifferent strand: $transcript_diff\n";
	print "\tTotal number UTR, CDS overlap: $counter_transcript_CDS_overlap\n";
	print "\t\tsame strand: $transcript_CDS_same\n\t\tdifferent strand: $transcript_CDS_diff\n";
	print "\tTotal possible interdigitating genes: $counter_interdigitates\n\t\tsame strand: $interdigitates_same\n\t\tdifferent strand: $interdigitates_diff\n\n";
}

#####################################
#			 SUBROUTINES 			#
#####################################
# example: &file_reader_alt ($gff_subfile, $dir)
# stores file info into hashes
sub file_reader_alt {
	my $subfile = shift;
	#empty hashes b4 iterating through a new file
	%gene_to_start = (); %gene_to_stop = (); %gene_to_strand = (); %gene_to_CDSs = ();
	%CDS_to_gene = (); %CDS_to_start = (); %CDS_to_stop = (); %CDS_to_exons = (); %CDS_to_introns = (); %CDS_to_strand = ();
	%exon_to_start = (); %exon_to_stop = ();
	%intron_to_start = (); %intron_to_stop = ();
	%transcript_to_start = (); %transcript_to_stop = (); %transcript_to_strand = (); %CDS_to_transcript = ();

	open (IN, "<$subfile") or die "file_reader subroutine can not open file name $subfile\n";
	while (my $line = <IN>){
		next if ($line =~ /^#/);		
		chomp $line;
		my ($chromosome, $feature, $featureid, $start, $stop, $strand, $comment) = split (/\t/, $line);
		#storing neccessary data into hashes for CDS
		if ($feature eq 'CDS') {
			my ($gid, $status) = split (/ ; /, $comment);
			$CDS_to_gene{$featureid} = $gid;					
			$CDS_to_start{$featureid} = $start;
			$CDS_to_stop{$featureid} = $stop;
			$CDS_to_strand{$featureid} = $strand;
			$gene_to_CDSs{$gid} .= "$featureid:";
		}
		#storing neccessary data into hash for exon
		if ($feature eq 'exon') {
			$exon_to_start{"$start-$stop"} = $start; 
			$exon_to_stop{"$start-$stop"} = $stop;
			$CDS_to_exons{$featureid} .= "$start,$stop:";
		}
		#storing neccessary data into hash for intron
		if ($feature eq 'intron') {
			$intron_to_start{"$start-$stop"} = $start; 
			$intron_to_stop{"$start-$stop"} = $stop;
			$CDS_to_introns{$featureid} .= "$start,$stop:";
		}
		#storing neccessary data into hashes for gene
		if ($feature eq 'gene') {
			$counter_gene++;
			$gene_to_start{$featureid} = $start;
			$gene_to_stop{$featureid} = $stop;
			$gene_to_strand{$featureid} = $strand;
		}
		#store neccessary data into hashes for gene
		if ($feature eq 'transcript') {
			$transcript_to_start{$featureid} = $start;
			$transcript_to_stop{$featureid} = $stop;
			$transcript_to_strand{$featureid} = $strand;
			$CDS_to_transcript{$comment} .= "$featureid:"; 
		}
	}
}

#check overlap in CDS
sub CDS_check {
	my ($geneA, $geneB) = @_; 
	#create key to look up what's returned in $a
	my %overlap_description = (
		0 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tno overlap between geneA and geneB",
		1 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\t3'end of geneA overlaps 5'end of geneB",
		2 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\t5'end of geneA overlaps 3' end of geneB",
		3 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tgeneB contains geneA",
		4 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tgeneA contains geneB",
		5 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tgeneA and geneB are identical",	
		6 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tgeneA contained in intron of geneB", #CDS/Intron
		7 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tgeneB contained in intron of geneA", #Intron/CDS
		8 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\texons overlap", 
		9 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tgeneA and geneB possibly interdigitates", #Gene/Gene
		10 => "geneA\tcdsA\tgeneB\tcdsB\tstrandA\/strandB\tno exon overlap, geneA and geneB interdigitates", 
	);
	#obtain list of CDS for gene
	my @CDS1 = split (/:/,$gene_to_CDSs{$geneA});
	my @CDS2 = split (/:/,$gene_to_CDSs{$geneB});
	
	#loop through CDS to find which pair overlaps
	for (my $k = 0; $k < @CDS1; $k++) {
		LABEL: for (my $l = 0; $l < @CDS2; $l++) {
			
			#test for overlapping CDSs
			my ($a,$b) = &overlap_check($CDS_to_start{$CDS1[$k]}, $CDS_to_stop{$CDS1[$k]}, $CDS_to_start{$CDS2[$l]}, $CDS_to_stop{$CDS2[$l]});
			
			if ($a == 0) { #cds doesn't overlap, overlap in UTRs? 
				# place cds in subroutine to check for overlap in UTR
				my ($types, $a1, $b1, $transcript1, $transcript2) = &transcript_overlap_check ($CDS1[$k], $CDS2[$l]);
				next LABEL unless ($a1 && $a1 != 0);
				#place the outputs into subroutine to generate the proper description from %overlap_description
				my $description_return = &description ($overlap_description{$a1}, $geneA, $geneB, $transcript1, $transcript2);
				return "$types\t$b1\t$description_return";
			}
			#test to see if CDS is in intron of 2nd gene
			my ($result, $comment) = &gene_in_intron_check ($CDS1[$k], $CDS2[$l]);	
			#move to next cds if neither has introns, has to have introns for interdigitation
			next LABEL if ($result eq 'next');
			
			#at this point, the CDSs overlaps, checks to see if CDS is contained wholly in an intron
			if ($result == 6 or $result == 7) {
				my $types;
				$types = "CDS/Intron" if $result == 6;
				$types = "Intron/CDS" if $result == 7;
				#place the outputs into subroutine to generate the proper description from %overlap_description
				my $description_return2 = &description ($overlap_description{$result}, $geneA, $geneB, $CDS1[$k], $CDS2[$l]);
				return "$types\t$comment\t$description_return2";
			}
			
			#if geneA isn't in intron of geneB or vice versa, then test for interdigitation
			elsif ($result != 6 or $result != 7) { 
				next LABEL if (!exists $CDS_to_introns{$CDS1[$k]} or !exists $CDS_to_introns{$CDS2[$l]});
				my ($result1) = &interdigitation_check ($CDS1[$k], $CDS2[$l]);
				if ($result1 == 9) {
					#place the outputs into subroutine to generate the proper description from %overlap_description
					my $description_return3 = &description ($overlap_description{$result1}, $geneA, $geneB, $CDS1[$k], $CDS2[$l]);		
					return "Gene/Gene\t\.\t$description_return3";		
				}
			}
		}
	}
	return "null";
}

#checks nonoverlaping CDS to see if transcript has overlap
sub transcript_overlap_check {

	my ($cds1, $cds2) = @_;
	
	# in the case one CDS does not have transcript records:
	# 1. CDS1 with UTR of CDS2
	if (exists $CDS_to_transcript{$cds2} and !exists $CDS_to_transcript{$cds1}) {
		my @transcript2 = split (/:/, $CDS_to_transcript{$cds2});
		foreach my $transcript2 (@transcript2) {
			my ($a,$b) = &overlap_check ($CDS_to_start{$cds1}, $CDS_to_stop{$cds1}, $transcript_to_start{$transcript2}, $transcript_to_stop{$transcript2});
			next if ($a == 0);
			$counter_transcript_CDS_overlap ++;
			$transcript_CDS_same++ if ($transcript_to_strand{$transcript2} eq $CDS_to_strand{$cds1});
			$transcript_CDS_diff++ if ($transcript_to_strand{$transcript2} ne $CDS_to_strand{$cds1});
			return ("3\'CDS/5\'UTR", $a, $b, $cds1, $transcript2) if ($a == 1);
			return ("5\'CDS/3\'UTR", $a, $b, $cds1, $transcript2) if ($a == 2);
			return ("CDS/UTR", $a, $b, $cds1, $transcript2)
		}
	}
	
	#2. CDS2 with UTR of CDS1
	elsif (exists $CDS_to_transcript{$cds1} and !exists $CDS_to_transcript{$cds2}) {
		my @transcript1 = split (/:/, $CDS_to_transcript{$cds1});
		foreach my $transcript1 (@transcript1) {
			my ($a,$b) = &overlap_check ($transcript_to_start{$transcript1}, $transcript_to_stop{$transcript1}, $CDS_to_start{$cds2}, $CDS_to_stop{$cds2});
			next if ($a == 0);
			$counter_transcript_CDS_overlap ++;
			$transcript_CDS_same++ if ($CDS_to_strand{$cds2} eq $transcript_to_strand{$transcript1});
			$transcript_CDS_diff++ if ($CDS_to_strand{$cds2} ne $transcript_to_strand{$transcript1});
			return ("3\'UTR/5\'CDS", $a, $b, $transcript1, $cds2) if ($a == 1);
			return ("5\'UTR/3\'CDS", $a, $b, $transcript1, $cds2) if ($a == 2);
			return ("UTR/CDS", $a, $b, $transcript1, $cds2)
		}
	}
	
	#both CDS has transcripts
	elsif (exists $CDS_to_transcript{$cds1} and exists $CDS_to_transcript{$cds2}) {

		my @transcript1 = split (/:/, $CDS_to_transcript{$cds1});
		my @transcript2 = split (/:/, $CDS_to_transcript{$cds2});
		
		foreach my $transcript1 (@transcript1){

			#compare cds2 with transcripts of cds1
			my ($code, $nucleotide) = &overlap_check ($transcript_to_start{$transcript1}, $transcript_to_stop{$transcript1}, $CDS_to_start{$cds2}, $CDS_to_stop{$cds2});
			if ($code != 0) {
				$counter_transcript_CDS_overlap ++;
				$transcript_CDS_same++ if ($CDS_to_strand{$cds2} eq $transcript_to_strand{$transcript1});
				$transcript_CDS_diff++ if ($CDS_to_strand{$cds2} ne $transcript_to_strand{$transcript1});
			}
			return ("3\'UTR/5\'CDS", $code, $nucleotide, $transcript1, $cds2) if ($code == 1);
			return ("5\'UTR/3\'CDS", $code, $nucleotide, $transcript1, $cds2) if ($code == 2);			
			return ("UTR/CDS", $code, $nucleotide, $transcript1, $cds2) if ($code != 0);

			foreach my $transcript2 (@transcript2){
			
				#compare cds1 with transcripts of cds2
				my ($a1,$b2) = &overlap_check ($CDS_to_start{$cds1}, $CDS_to_stop{$cds1}, $transcript_to_start{$transcript2}, $transcript_to_stop{$transcript2});
				if ($a1 != 0) {
					$counter_transcript_CDS_overlap++;
					$transcript_CDS_same++ if ($CDS_to_strand{$cds1} eq $transcript_to_strand{$transcript2});
					$transcript_CDS_diff++ if ($CDS_to_strand{$cds1} ne $transcript_to_strand{$transcript2});
				}
				return ("3\'CDS/5\'UTR", $a1, $b2, $cds1, $transcript2) if ($a1 == 1);
				return ("5\'CDS/3\'UTR", $a1, $b2, $cds1, $transcript2)	if ($a1 == 2);			
				return ("CDS/UTR", $a1, $b2, $cds1, $transcript2) if ($a1 != 0);
				
				#if CDS and UTR doesn't overlap, check for overlapping UTRs
				my ($a,$b) = &overlap_check ($transcript_to_start{$transcript1}, $transcript_to_stop{$transcript1}, $transcript_to_start{$transcript2}, $transcript_to_stop{$transcript2});
				next if ($a == 0);
				$counter_transcript_overlap ++;
				$transcript_same++ if ($transcript_to_strand{$transcript2} eq $transcript_to_strand{$transcript1});
				$transcript_diff++ if ($transcript_to_strand{$transcript2} ne $transcript_to_strand{$transcript1});
				return ("3\'UTR/5\'UTR", $a, $b, $transcript1, $transcript2) if ($a == 1);
				return ("5\'UTR/5\'UTR", $a, $b, $transcript1, $transcript2) if ($a == 2);				
				return ("UTR/UTR", $a, $b, $transcript1, $transcript2) if ($a != 2 or $a != 1);
			}
		}
	}
	return 0;
}

#gets proper description for RNA gene overlap
sub description1 {
	my ($string, $geneA, $geneB) = @_;
	if ($web) {
		my $url = $base_url . "$geneA" . ";class=";
		$string =~ s/strandA/$gene_to_strand{$geneA}/;
		$string =~ s/strandB/$gene_to_strand{$geneB}/;
		$string =~ s/geneA\tcdsA\tgeneB\tcdsB\t//;
		$string =~ s/geneA/$geneA/;
		$string =~ s/geneB/$geneB/;	
		$string .= "\t$url";
	}
	$string =~ s/geneA/$geneA/g if (!$web);
	$string =~ s/geneB/$geneB/g if (!$web);	
	$string =~ s/cdsA/""/;
	$string =~ s/cdsB/""/;
	$string =~ s/strandA/$gene_to_strand{$geneA}/;
	$string =~ s/strandB/$gene_to_strand{$geneB}/;
	return $string;
}				

#gets proper description for protein coding gene overlap
sub description {
	my ($string, $geneA, $geneB, $cds1, $cds2) = @_;
	if ($web) {
		my $url = $base_url . "$geneA" . ";class=";
		$string =~ s/strandA/$gene_to_strand{$geneA}/;
		$string =~ s/strandB/$gene_to_strand{$geneB}/;
		$string =~ s/geneA\tcdsA\tgeneB\tcdsB\t//;
		$string =~ s/geneA/$geneA/;
		$string =~ s/geneB/$geneB/;	
		$string .= "\t$url";
	}			
	$string =~ s/geneA/$geneA/g if (!$web);
	$string =~ s/cdsA/$cds1/ if (!$web);
	$string =~ s/cdsB/$cds2/ if (!$web);
	$string =~ s/geneB/$geneB/g if (!$web);
	$string =~ s/strandA/$gene_to_strand{$geneA}/;
	$string =~ s/strandB/$gene_to_strand{$geneB}/;
	return $string; # returns WBGene00001	CDS1	WBGene00002	CDS2	+/-	gene in intron of gene
}

#put in genes and their CDSs, checks if gene is located in intron of the other gene
sub gene_in_intron_check {
	my ($cdsA, $cdsB) = @_;
	
	#need to make sure at least one of the CDS has a list of introns
	#else exit the subroutine and go to next cds (print out the CDS overlap type in all option?)
	return ("next", "next") if (!exists $CDS_to_introns{$cdsA} and !exists $CDS_to_introns{$cdsB});

	#grab all introns belonging to CDS containing introns
	my @intronsA = split(/:/, $CDS_to_introns{$cdsA}) if (exists $CDS_to_introns{$cdsA});
	my @intronsB = split(/:/, $CDS_to_introns{$cdsB}) if (exists $CDS_to_introns{$cdsB});
	
	#if cdsB doesn't have intron... 
	if (!exists $CDS_to_introns{$cdsB}) {
		#loop through cdsA to test if whole cdsB is in any intron of cdsA
		foreach my $intronA (@intronsA){
			$intronA =~ m/(\d+),(\d+)/;
			my ($startA, $stopA) = ($1,$2);

			# we're only interested in cdsB in intron of cdsA:
			# overlap_check option 3: E(3/4) contains E(1/2)			
			# more specifically, in this case, option 3 points to overlap_description 7 => "geneB contained in intron of geneA\t"
			my ($a, $b) = &overlap_check($CDS_to_start{$cdsB}, $CDS_to_stop{$cdsB}, $startA, $stopA);
			 if ($a == 3) {
				$counter_gene_in_intron++;
				$intron_same++ if ($CDS_to_strand{$cdsA} eq $CDS_to_strand{$cdsB});
				$intron_diff++ if ($CDS_to_strand{$cdsA} ne $CDS_to_strand{$cdsB});
			}
			return (7, $b) if ($a == 3);
		}
	}

	#if cdsA doesn't have intron... 
	elsif (!exists $CDS_to_introns{$cdsA}) {
		#loop through cdsA to test if whole cdsB is in any intron cdsA
		foreach my $intronsB (@intronsB){
			$intronsB =~ m/(\d+),(\d+)/;
			my ($startB, $stopB) = ($1,$2);
			# we're only interested in cdsA in intron of cdsB:
			# overlap_option 3: E(3/4) contains E(1/2)
			# more specifically, in this case, option 3 points to overlap_description 6 => "geneA contained in intron of geneB\t",
			my ($a1, $b1) = &overlap_check($CDS_to_start{$cdsA}, $CDS_to_stop{$cdsA}, $startB, $stopB);
			if ($a1 == 3) {
				$counter_gene_in_intron++;
				$intron_same++ if ($CDS_to_strand{$cdsA} eq $CDS_to_strand{$cdsB});
				$intron_diff++ if ($CDS_to_strand{$cdsA} ne $CDS_to_strand{$cdsB});
			}
			return (6, $b1) if ($a1 == 3);
		}
	}
	
	#if both cds has introns, find to see if gene/CDS is located inside the intron
	else {
		foreach my $introns_A (@intronsA) {
			
			$introns_A =~ m/(\d+),(\d+)/;
			my ($startA1, $stopA1) = ($1,$2);
			
			foreach my $introns_B (@intronsB) {
				
				$introns_B =~ m/(\d+),(\d+)/;
				my ($startB1, $stopB1) = ($1,$2);
				
				#compare the genes to the introns
				my ($a_1, $b_1) = &overlap_check($CDS_to_start{$cdsA}, $CDS_to_stop{$cdsA}, $startB1, $stopB1);
				my ($a_2, $b_2) = &overlap_check($CDS_to_start{$cdsB}, $CDS_to_stop{$cdsB}, $startA1, $stopA1);
			
				#report if genes are located in intron (3. E(3/4) contains E(1/2))
				if ($a_1 == 3 or $a_2 == 3) {
					$counter_gene_in_intron++;
					$intron_same++ if ($CDS_to_strand{$cdsA} eq $CDS_to_strand{$cdsB});
					$intron_diff++ if ($CDS_to_strand{$cdsA} ne $CDS_to_strand{$cdsB});
				}
				#returns a number code to be look up and the amount of nucleotide overlap
				return (6, abs($CDS_to_start{$cdsA}-$CDS_to_stop{$cdsA})) if ($a_1 == 3);
				return (7, abs($CDS_to_start{$cdsB}-$CDS_to_stop{$cdsB})) if ($a_2 == 3);		
			}
		}
	}
}

#put in two genes and their CDSs, get out 2 lists of exons, then test for intergitation
sub interdigitation_check{	
	my ($cdsA,$cdsB) = @_;
	
	#grab all exons belonging to each CDS
	my @exonsA = split(/:/, $CDS_to_exons{$cdsA}); 
	my @exonsB = split(/:/, $CDS_to_exons{$cdsB});
	
	#looping through the exons from the exon lists to check for overlap
	foreach my $exonsA (@exonsA){
	
		$exonsA =~ m/(\d+),(\d+)/;
		my ($startA,$stopA) = ($1,$2);
		
		foreach my $exonsB (@exonsB){
	
			$exonsB =~ m/(\d+),(\d+)/;
			my ($startB,$stopB) = ($1,$2);
			
			#check exons against exons
			my ($overlap_value_return, $b3) = &overlap_check($startA,$stopA,$startB,$stopB);
		
			return (8) if ($overlap_value_return != 0); #exons overlap
			if ($overlap_value_return ==0) {
				$counter_interdigitates++;
				$interdigitates_same++ if ($CDS_to_strand{$cdsA} eq $CDS_to_strand{$cdsB});
				$interdigitates_diff++ if ($CDS_to_strand{$cdsA} ne $CDS_to_strand{$cdsB});
			}
			return (9) if ($overlap_value_return == 0); #exons does not overlap, possibly interdigitates
		}
	}
}

#compares 2 sets of start/stop coordinates
sub overlap_check {

	# start/stop of 1st set
	my ($xstart, $xstop, $ystart, $ystop) = @_;

	###### 0. eliminate all non-matching genes (A and B does not match)######
	if	(($ystop < $xstart) or ($xstop < $ystart)) {return (0, 'NA');}
				
	###### 1. 3'end of E(1/2) overlaps 5' end of E(3/4) ######
	elsif 	(
				(($xstart < $ystart) and ($ystart<$xstop) and ($xstop<$ystop)) or
				(($xstart < $ystart) and ($ystart<=$xstop) and ($xstop<$ystop))
			)
			{return (1, abs($xstop - $ystart));}
				
	###### 2. 5'end of E(1/2) overlaps 3' end of E(3/4) ######
	elsif 	(
				(($ystart<$xstart) and ($xstart<$ystop) and ($ystop<$xstop)) or 
				(($ystart<$xstart) and ($xstart<=$ystop) and ($ystop<$xstop))
			)
			{return (2, abs($ystop - $xstop));}
				
	###### 3. E(3/4) contains E(1/2) ######
	elsif 	(
				(($ystart < $xstart) and ($xstop < $ystop)) or
				(($ystart <= $xstart) and ($xstop < $ystop)) or
				(($ystart < $xstart) and ($xstop <= $ystop))
			)
			{return (3, abs($xstop - $xstart));}
			
	###### 4. E(1/2) contains E(3/4) ######		
	elsif	(
				(($xstart < $ystart) and ($ystop < $xstop)) or
				(($xstart <= $ystart) and ($ystop < $xstop)) or
				(($xstart < $ystart) and ($ystop <= $xstop))
			)
			{return (4, abs($ystop - $ystart));}
				
	###### 5. E(1/2) identical to E(3/4) ######
	elsif	(($xstart == $ystart) and ($xstop == $ystop))
			{return (5, $xstop - $xstart);}
			
	###### 6. makes sure I didn't miss any other possibilities ######
	else	{die "logic FAILS!: ($xstart,$xstop,$ystart,$ystop)\n";}
}