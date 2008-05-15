#!/usr/bin/perl
#
# make_bin_and_lib_dirs.pl
#
# A script to make symbolic links of perl scripts (*.pl) and perl modules (*.pm) into a new bin and lib directory
#
# by Keith Bradnam, May 2008
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use Getopt::Long;
use File::Glob ':glob';


########################
# Command line options
########################

my $source; # directory containing original files
my $target; # directory where sym links will go
my $pl;     # create sym links for perl scripts (*.pl)
my $pm;     # create sym links for perl modules (*.pm)
 
GetOptions ("source=s"  => \$source,
			"target=s"  => \$target,
			 "pl"       => \$pl,
			 "pm"       => \$pm);

# check command line options
die "must use -source option to specify a directory that contains some perl scripts or modules\n" if (!$source);
die "must specify -target option for location of directory where symbolic links will be created\n" if (!$target);
die "-source option has to specify full path to source directory e.g. ~keith/perl or /Users/keith/perl\n" if ($source !~ m/^(\/|~)/);
die "must specify -pl to make links for perl scripts, -pm to make links for perl modules or both -pl and -pm\n" if (!$pl && !$pm);


# read source directory and get list of perl scripts or modules
# using glob module to get paths
my @files;
@files = bsd_glob("$source/*.pl")    if ($pl  && !$pm);
@files = bsd_glob("$source/*.pm")    if (!$pl && $pm);
@files = bsd_glob("$source/*.p[lm]") if ($pl  && $pm);


# make symbolic links
foreach my $file (@files){
	system("ln -s $file $target\n") && print "\nWarning: Couldn't create symbolic link for $file\n";
}

exit(0);




