#!/usr/bin/perl
#
# tandemify.pl
#
# A script to make a tandem version of each sequence within a FASTA file
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict; use warnings;
use Keith;


# a short script as all the work is done by the function in Keith.pm
die "Usage: $0 <fasta file> <optional output file name>\n" unless (@ARGV == 1 or @ARGV == 2);

my ($input, $output) = @ARGV; 
Keith::tandemify_sequence($input, $output);

exit(0);

