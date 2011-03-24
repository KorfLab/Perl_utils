#!/usr/bin/perl
use strict; use warnings;
use FAlite;
use Getopt::Std;
use vars qw($opt_d $opt_g $opt_h $opt_s $opt_w $opt_W $opt_x $opt_i $opt_G $opt_T $opt_m $opt_n $opt_q $opt_r);
getopts('dgh:s:w:W:x:i:GT:m:n:q:r:');

my $MINSCORE  = 20;
my $MAXCONTIG = 100000;
my $WORDSIZE  = 13;
my $WINK      = 0;
my $HSPMAX    = 0;
my $PERCENT   = 75;
my $GLOBALT   = 2;
my $MATCH     = 1;
my $MISMATCH  = -1;
my $GAPINS    = 2;
my $GAPEXT    = 2;

die "
QSTACK - Multiplex BLASTN (variable read length)
usage: qstack.pl [options] <blastn db> <fasta file>
options:
  -d           dust sequence         [default off]
  -g           allow gaps            [default no gaps]
  -h  <int>    hsp max               [default $HSPMAX]
  -s  <int>    min score             [default $MINSCORE, +1,-1,-2 method]
  -i  <float>  min percent           [default $PERCENT]
  -w  <int>    wordsize              [default $WORDSIZE]
  -W  <int>    wink                  [default unspecified]
  -x  <int>    multiplex length      [default $MAXCONTIG]
  -G           global-ish align      [default off]
  -T  <int>    global threshold      [default $GLOBALT]
  -m  <int>    match score           [default $MATCH]
  -n  <int>    mismatch score        [default $MISMATCH]
  -q  <int>    gap insertion penalty [default $GAPINS]
  -r  <int>    gap extension penalty [default $GAPEXT]

" unless @ARGV == 2;

my ($DB, $QUERY) = @ARGV;
my $DUST   = $opt_d ? 'filter=dust' : '';
$GAPINS    = $opt_q if $opt_q;
$GAPEXT    = $opt_r if $opt_r;
my $GAPS   = $opt_g ? "Q=$GAPINS R=$GAPEXT" : 'nogap';
$MINSCORE  = $opt_s if $opt_s;
$WORDSIZE  = $opt_w if $opt_w;
$WINK      = $opt_W if $opt_W;
$MAXCONTIG = $opt_x if $opt_x;
$HSPMAX    = $opt_h if $opt_h;
$PERCENT   = $opt_i if $opt_i;
my $GLOBAL = $opt_G;
$GLOBALT   = $opt_T if $opt_T;
$MATCH     = $opt_m if $opt_m;
$MISMATCH  = $opt_n if $opt_n;

###############################################################################
# Setup
###############################################################################

my $tmp_dna = "tmp.qstack.$$.dna";
my $tmp_rep = "tmp.qstack.$$.blast";
my $p1 = "M=$MATCH N=$MISMATCH kap mformat=3 hspmax=0 B=2147483647 V=0"; # warnings? notes?
my $p2 = "$GAPS W=$WORDSIZE S=$MINSCORE S2=$MINSCORE";
$p2 .= " hspmax=$HSPMAX -warnings" if $HSPMAX;
$p2 .= " WINK=$WINK" if $WINK;
open(IN, $QUERY) or die;
my $fasta = new FAlite(\*IN);

###############################################################################
# Main Loop
###############################################################################

my $contig_length = 0;
my @contig;
my @map;
my %seq_length;

while (my $entry = $fasta->nextEntry) {

	# contig
	my ($def) = $entry->def =~ /^>(\S+)/;
	my $seq = $entry->seq;
	$seq_length{$def} = length($seq);
	push @contig, $seq;
	my $length = length($seq) +1;
	for (my $i = $contig_length+1; $i <= $contig_length + $length; $i++) {
		$map[$i] = {
			def => $def,
			off => $contig_length,
		}
	}
	$contig_length += $length;
	next unless $contig_length > $MAXCONTIG;
	
	# blast	
	run_blast();
	reformat_blast();
	
	# reset
	$contig_length = 0;
	@map = ();
	@contig = ();
}

if (@contig) {
	run_blast();
	reformat_blast();
}

sub run_blast {
	open(OUT, ">$tmp_dna") or die;
	print OUT ">qstack-sequence\n";
	foreach my $contig (@contig) {
		print OUT $contig, "\n", "-", "\n";
	}
	close OUT;
	system("blastn $DB $tmp_dna $p1 $p2 > $tmp_rep") == 0 or die;
}

sub reformat_blast {
	open(BLAST, $tmp_rep) or die;
	while (<BLAST>) {
		next unless /^\w/;
		my @f = split;
		my $off = $map[$f[17]]{off};
		my $def = $map[$f[17]]{def};
		if (not defined $off) {
			use DataBrowser; browse(\@map);exit
		}
		$f[17] -= $off;
		$f[18] -= $off;
		$f[0]   = $def;
		next if $f[10] < $PERCENT;
		if ($GLOBAL) {
			my $align_length = abs($f[17] - $f[18]) + 1;
			next if $seq_length{$def} - $align_length > $GLOBALT;
		}
		print join("\t", @f), "\n";
	}
	close BLAST;
}

###############################################################################
# Clean up
###############################################################################

END {
	unlink $tmp_dna if defined $tmp_dna and -e $tmp_dna;
	unlink $tmp_rep if defined $tmp_rep and -e $tmp_rep;
}
