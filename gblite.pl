#!/usr/bin/perl
###############################################################################
# Submitted to the Public Domain by author Ian Korf
# Version 2008-02-08
###############################################################################
use strict; use warnings; use sigtrap;
use FileHandle;
use IPC::Open2;
use Getopt::Std;
use vars qw($opt_p $opt_a $opt_m $opt_t $opt_T $opt_c $opt_e);
use threads; use Thread::Queue;
my $usage = "usage: gblite.pl [options] <command>
commands:
  build  <gb> <tax> <build>    process all files
    -p <int>                   parsing threads [default 2]
  search
    -a                         amino acids [default is nucleotides]
    -m <string>                nucleic acids with this moltype
    -t <taxid>                 only sequences from this taxid
    -T <taxid>                 omit sequences from this taxid
    -c                         just count, do not retrieve sequences
    -e                         echo SQL commands\n";
die $usage unless @ARGV;

# options
getopts('p:am:t:T:ce');
my $THREADS    = $opt_p ? $opt_p : 2;
my $COUNT_ONLY = $opt_c;
my $ECHO_SQL   = $opt_e;
my $PROTEIN    = $opt_a;
my $MOLTYPE    = $opt_m;
my $TAXNODE    = $opt_t;
my $TAXNOT     = $opt_T;

# environment
my $ROOT = $ENV{GBLITE};
die "GBLITE must point to directory\n" unless defined $ROOT and -d $ROOT;
my $DB = "$ROOT/gblite.db";

# command processing
my $CMD = shift @ARGV;
if    ($CMD eq 'build')  {build(@ARGV)}
elsif ($CMD eq 'parse')  {parse_gb(@ARGV)}
elsif ($CMD eq 'search') {search()}
else                     {die "unrecognized command: $CMD\n$usage"}
exit(0);

###############################################################################
# Subroutines
###############################################################################

sub load_table {
	my ($file, $table) = @_;
	print "importing $file into $table table\n";
	open(SQL, "| sqlite3 $DB") or die;
	print SQL '.separator \t', "\n";
	my ($tmp) = "/tmp/gblite.$$.tmp";
	system("gunzip -c $file > $tmp") == 0 or die;
	print SQL ".import $tmp $table\n";
	close SQL;
	unlink $tmp;
}

sub worker {
	my ($q, $BUILD) = @_;
	my $tid = threads->tid;
	while ($q->pending) {
		my $file = $q->dequeue;
		print "processing file $file in thread $tid\n";
		system("gblite.pl parse $file $BUILD") == 0 or die;
	}
}

sub build {
	my ($GB, $TAX, $BUILD) = @_;
	die $usage unless -d $GB and -d $TAX and -d $BUILD;
		
	if (-e $DB) {
		print "Previous database exists. [return] to rebuild\n"; <STDIN>;
		unlink $DB;
	}
	
	open(SQL, "| sqlite3 $DB") or die;
	print SQL "CREATE TABLE GenBank (accession TEXT, taxid INT, moltype TEXT, source TEXT, PRIMARY KEY (accession));\n";
	print SQL "CREATE TABLE Taxonomy (taxid INT, name TEXT, parent INT, idxL INT, idxR INT, PRIMARY KEY (taxid));\n";
	close SQL;
	
	parse_tax("$TAX/names.dmp", "$TAX/nodes.dmp", $BUILD);
	my @gbfile = `ls $GB/*.seq.gz`; chomp @gbfile;
	my $Q = new Thread::Queue;
	$Q->enqueue(@gbfile);
	my @worker;
	for (my $i = 0; $i < $THREADS; $i++) {$worker[$i] = threads->create(\&worker, $Q, $BUILD)}
	for (my $i = 0; $i < $THREADS; $i++) {$worker[$i]->join}
	
	load_table("$BUILD/taxonomy.gz", "Taxonomy");
	my @table = `ls $BUILD/*.tbl.gz`; chomp @table;
	foreach my $table (@table) {load_table($table, 'GenBank')}
	system("xdformat -n -I -C N -o $ROOT/gblite $BUILD/*.nt.gz") == 0 or die;
	system("xdformat -p -I -C X -o $ROOT/gblite $BUILD/*.aa.gz") == 0 or die;
	print "\n*** gblite build complete ****\n";
}

sub nextEntry {
	my ($FH) = @_;
	my %e;
	while (<$FH>) {last if /^LOCUS/}
	return 0 if not defined $_;
	my @f = split;
	($e{acc}, $e{type}) = ($f[1], $f[4]);
	while (<$FH>) {last if /^DEFINITION/}
	chomp;
	$e{def} = substr($_, 12);
	while (<IN>) {
		if (/^\s+(\S.+)/) {$e{def} .= $1}
		else {last}
	}
	while (<$FH>) {last if /^FEATURES/}
	@f = ();
	while (<IN>) {
		last if /^ORIGIN|^CONTIG/;
		push @f, $_;
	}
			
	# taxid
	foreach (@f) {
		if (/^\s+\/db_xref="taxon:(\d+)"/) {$e{taxid} = $1; last}
	}
	$e{taxid} = -1 if not defined $e{taxid};

	# protein
	my (%cds, $pid, $pep, $product);
	my $translating = 0;
	for(@f) {
		if (/^\s+\/product="(.+)/) {
			$product = $1;
			$product =~ s/"//g;
		} elsif (/^\s+\/protein_id="(\w+)/) {
			$pid = $1;
		} elsif (/^\s+\/translation="(\w+)/) {
			$cds{$pid}{pep} = $1;
			$cds{$pid}{product} = $product;
			$translating = 1;
		} elsif (/^\s+\// or /^\s{5}\S+/) {
			$translating = 0;
		} elsif (/^\s+(\w+)/ and $translating) {
			$cds{$pid}{pep} .= $1
		}
	}
	$e{cds} = \%cds;
	
	# dna
	$e{dna} = "";
	while (<$FH>) {
		last if /^\/\//;
		my @s = split;
		shift @s;
		$e{dna} .= join("", @s);
	}
	
	return \%e;
}

sub parse_gb {
	my ($file, $BUILD) = @_;
	my ($name) = $file =~ /gb(\w+)\.seq\.gz/;
	my $NTFILE = "$BUILD/$name.nt.gz";
	my $AAFILE = "$BUILD/$name.aa.gz";
	my $TABLE  = "$BUILD/$name.tbl.gz";
	return if -e $NTFILE and -e $AAFILE and -e $TABLE;
	open(NOUT, "| gzip > $NTFILE") or die;
	open(POUT, "| gzip > $AAFILE") or die;
	open(TOUT, "| gzip > $TABLE") or die;
	open(IN, "gunzip -c $file |") or die;
	while (my $e = nextEntry(\*IN)) {
		my %cds = %{$e->{cds}};
		foreach my $id (keys %cds) {
			my $desc = defined $cds{$id}{product} ? $cds{$id}{product} : "";
			print POUT ">", $id, " ", $desc, " taxid:", $e->{taxid},  " source:", $e->{acc}, "\n", $cds{$id}{pep}, "\n";
			print TOUT $id, "\t", $e->{taxid}, "\taa\t", $e->{acc}, "\n";
		}
		if ($e->{dna}) {
			print NOUT ">", $e->{acc}, " ", $e->{def}, " taxid:", $e->{taxid}, "\n", $e->{dna}, "\n";
			print TOUT $e->{acc}, "\t", $e->{taxid}, "\t", $e->{type}, "\t\n";
		}
		
	}
	close IN; close NOUT; close POUT; close TOUT;
}

sub createTree {
	my ($P, $C, $I, $int, $parent, $tree) = @_;
	foreach my $child (keys %{$C->{$parent}}) {
		$tree->{$parent}{$child} = {};
		$I->{$child}{LEFT} = ++$$int;
		createTree($P, $C, $I, $int, $child, $tree->{$parent});
		$I->{$child}{RIGHT} = ++$$int;
	}
	return $tree;
}

sub parse_tax {
	my ($fnames, $fnodes, $BUILD) = @_;
	print "processing $fnames and $fnodes\n";
	return if -e "$BUILD/taxonomy.gz";
	my %SN;
	open(IN, $fnames) or die "$fnames not found";
	while (<IN>) {
		next unless /scientific name/;
		my ($taxid, $name) = split(/\|/, $_);
		$name =~ s/^\s+//g;
		$name =~ s/\s+$//g;
		$taxid =~ s/\s//g;
		$SN{$taxid} = $name;
	}
	close IN;
	
	my (%Parent, %Child);
	open(IN, $fnodes) or die;
	while (<IN>) {
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
	
	my %Index;
	my $int = 1;
	my $root = 1; 
	my $tree = createTree(\%Parent, \%Child, \%Index, \$int, $root, {});
	$Index{1}{LEFT} = 0;
	$Index{1}{RIGHT} = ++$int;
	$Parent{1} = 0;
		
	open(OUT, "| gzip > $BUILD/taxonomy.gz") or die;
	foreach my $taxid (keys %SN) {
		print OUT join("\t", $taxid, $SN{$taxid}, $Parent{$taxid}, $Index{$taxid}{LEFT}, $Index{$taxid}{RIGHT}), "\n";
	}	
	close OUT;

}

sub tax_indexes {
	my ($node) = @_;
	my ($R, $W) = (new FileHandle, new FileHandle);
	my $pid = open2($R, $W, "sqlite3 $DB") or die;
	if ($node =~ /^\d+$/) {print $W "SELECT idxL, idxR FROM taxonomy WHERE taxid = $node;\n"}
	else {print $W "SELECT idxL, idxR from taxonomy WHERE name = \"$node\";\n"}
	my $line = <$R>;
	my ($idxL, $idxR) = $line =~ /^(\d+)\|(\d+)/;
	close $R; close $W;
	return $idxL, $idxR;
}

sub search {
	my $sql = "SELECT accession FROM GenBank ";
	if ($TAXNODE or $TAXNOT) {
		$sql .= "JOIN Taxonomy ON GenBank.taxid = Taxonomy.taxid WHERE ";
		if ($TAXNODE) {
			my ($idxL, $idxR) = tax_indexes($TAXNODE);
		 	$sql .= "(Taxonomy.idxL >= $idxL AND Taxonomy.idxR <= $idxR) AND ";
		 }
		 if ($TAXNOT) {
		 	my ($idxL, $idxR) = tax_indexes($TAXNOT);
		 	$sql .= "NOT (Taxonomy.idxL >= $idxL AND Taxonomy.idxR <= $idxR) AND ";
		 }
	} else {$sql .= "WHERE "}
	my @where;
	if    ($PROTEIN) {push @where, "GenBank.moltype = \"aa\""}
	elsif ($MOLTYPE) {push @where, "GenBank.moltype = \"$MOLTYPE\""}
	else             {push @where, "GenBank.moltype != \"aa\""}	
	$sql .= join(" AND ", @where) . ";\n";
	if ($COUNT_ONLY) {$sql =~ s/accession/count(accession)/}
	my $tmpfile = "/tmp/gblite.$$.txt";
	open(SQL, "| sqlite3 $DB") or die;
	print SQL ".output $tmpfile\n" unless $COUNT_ONLY;
	print SQL $sql, "\n";
	print STDERR $sql, "\n" if $ECHO_SQL;
	close SQL;
	return if $COUNT_ONLY;
	
	my $np = $PROTEIN ? "-p" : "-n";
	system("xdget $np -f $ROOT/gblite $tmpfile");
	unlink $tmpfile;
}
