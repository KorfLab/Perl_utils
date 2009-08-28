#!/usr/bin/perl
# this script is to run GFFlite.pm to create subset files of the worm genome gff_files
use strict; use warnings; 
use Getopt::Long; 
use worm_genome_GFFlite;

#user must specify a directory to look up
die "usage: $0 <dir>" unless @ARGV;

########## 1. COMMAND LINE OPTIONS ##########
my $overwrite;	# creates new subsetfile that will overwrite the previous subsetfile
GetOptions (
	"overwrite!" => \$overwrite,
);

######## 2.a. ACCESSING DIRECTORY ##########
my $dir = $ARGV[0];
my @filtered_files = (); # example: subfile_file.gff
opendir (DIR, "$dir") or die "can not open directory $dir\n";
if ($overwrite) {
	while ($_ = readdir(DIR)) {
		if ($_ =~ /\S+\.gff/ and $_ !~ /subfile_(\S+\.gff)/) { 
			print "overwriting file $_...";
			push (@filtered_files, GFFlite::gff_reader_alt($dir , $_));
			print "finished overwriting\n" ;
		}
	}
	close DIR;
}
my @unfiltered_files; # example: file.gff
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
			my $dir_output;
			$dir_output = toygff::gff_reader_alt ($dir , $var);
			push (@filtered_files, $dir_output);
			print "\tfile $var_filtered created\n";
		}
	}
}
if (!$overwrite) {foreach my $file (@filtered_files) {print "$file\n";}}
print "GFFlite_reader complete\n";