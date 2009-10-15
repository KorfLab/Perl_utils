package IK;
use strict; use warnings;

my %AA = (
	'AAA' => 'K',	'AAC' => 'N',	'AAG' => 'K',	'AAT' => 'N',
	'AAR' => 'K',	'AAY' => 'N',	'ACA' => 'T',	'ACC' => 'T',
	'ACG' => 'T',	'ACT' => 'T',	'ACR' => 'T',	'ACY' => 'T',
	'ACK' => 'T',	'ACM' => 'T',	'ACW' => 'T',	'ACS' => 'T',
	'ACB' => 'T',	'ACD' => 'T',	'ACH' => 'T',	'ACV' => 'T',
	'ACN' => 'T',	'AGA' => 'R',	'AGC' => 'S',	'AGG' => 'R',
	'AGT' => 'S',	'AGR' => 'R',	'AGY' => 'S',	'ATA' => 'I',
	'ATC' => 'I',	'ATG' => 'M',	'ATT' => 'I',	'ATY' => 'I',
	'ATM' => 'I',	'ATW' => 'I',	'ATH' => 'I',	'CAA' => 'Q',
	'CAC' => 'H',	'CAG' => 'Q',	'CAT' => 'H',	'CAR' => 'Q',
	'CAY' => 'H',	'CCA' => 'P',	'CCC' => 'P',	'CCG' => 'P',
	'CCT' => 'P',	'CCR' => 'P',	'CCY' => 'P',	'CCK' => 'P',
	'CCM' => 'P',	'CCW' => 'P',	'CCS' => 'P',	'CCB' => 'P',
	'CCD' => 'P',	'CCH' => 'P',	'CCV' => 'P',	'CCN' => 'P',
	'CGA' => 'R',	'CGC' => 'R',	'CGG' => 'R',	'CGT' => 'R',
	'CGR' => 'R',	'CGY' => 'R',	'CGK' => 'R',	'CGM' => 'R',
	'CGW' => 'R',	'CGS' => 'R',	'CGB' => 'R',	'CGD' => 'R',
	'CGH' => 'R',	'CGV' => 'R',	'CGN' => 'R',	'CTA' => 'L',
	'CTC' => 'L',	'CTG' => 'L',	'CTT' => 'L',	'CTR' => 'L',
	'CTY' => 'L',	'CTK' => 'L',	'CTM' => 'L',	'CTW' => 'L',
	'CTS' => 'L',	'CTB' => 'L',	'CTD' => 'L',	'CTH' => 'L',
	'CTV' => 'L',	'CTN' => 'L',	'GAA' => 'E',	'GAC' => 'D',
	'GAG' => 'E',	'GAT' => 'D',	'GAR' => 'E',	'GAY' => 'D',
	'GCA' => 'A',	'GCC' => 'A',	'GCG' => 'A',	'GCT' => 'A',
	'GCR' => 'A',	'GCY' => 'A',	'GCK' => 'A',	'GCM' => 'A',
	'GCW' => 'A',	'GCS' => 'A',	'GCB' => 'A',	'GCD' => 'A',
	'GCH' => 'A',	'GCV' => 'A',	'GCN' => 'A',	'GGA' => 'G',
	'GGC' => 'G',	'GGG' => 'G',	'GGT' => 'G',	'GGR' => 'G',
	'GGY' => 'G',	'GGK' => 'G',	'GGM' => 'G',	'GGW' => 'G',
	'GGS' => 'G',	'GGB' => 'G',	'GGD' => 'G',	'GGH' => 'G',
	'GGV' => 'G',	'GGN' => 'G',	'GTA' => 'V',	'GTC' => 'V',
	'GTG' => 'V',	'GTT' => 'V',	'GTR' => 'V',	'GTY' => 'V',
	'GTK' => 'V',	'GTM' => 'V',	'GTW' => 'V',	'GTS' => 'V',
	'GTB' => 'V',	'GTD' => 'V',	'GTH' => 'V',	'GTV' => 'V',
	'GTN' => 'V',	'TAA' => '*',	'TAC' => 'Y',	'TAG' => '*',
	'TAT' => 'Y',	'TAR' => '*',	'TAY' => 'Y',	'TCA' => 'S',
	'TCC' => 'S',	'TCG' => 'S',	'TCT' => 'S',	'TCR' => 'S',
	'TCY' => 'S',	'TCK' => 'S',	'TCM' => 'S',	'TCW' => 'S',
	'TCS' => 'S',	'TCB' => 'S',	'TCD' => 'S',	'TCH' => 'S',
	'TCV' => 'S',	'TCN' => 'S',	'TGA' => '*',	'TGC' => 'C',
	'TGG' => 'W',	'TGT' => 'C',	'TGY' => 'C',	'TTA' => 'L',
	'TTC' => 'F',	'TTG' => 'L',	'TTT' => 'F',	'TTR' => 'L',
	'TTY' => 'F',	'TRA' => '*',	'YTA' => 'L',	'YTG' => 'L',
	'YTR' => 'L',	'MGA' => 'R',	'MGG' => 'R',	'MGR' => 'R',
);

sub complement {
	my ($seq) = @_;
	$seq =~ tr[ACGTRYMKWSBDHVacgtrymkwsbdhv]
	          [TGCAYRKMSWVHDBtgcayrkmwsvhdb];
	return $seq;
}

sub reverse_complement {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr[ACGTRYMKWSBDHVacgtrymkswbdhv]
	          [TGCAYRKMSWVHDBtgcayrkmswvhdb];
	return $seq;
}

sub kl_distance {
	my ($f1, $f2) = @_;
	
	my $D = 0;
	foreach my $k (keys %$f1) {
		$D += $f1->{$k} * log($f1->{$k} / $f2->{$k});
	}
	
	return $D;
}

sub count_to_freq {
	my ($hashref) = @_;
	
	my $total = 0;
	foreach my $k (keys %$hashref) {
		$total += $hashref->{$k};
	}
	
	my %freq;
	foreach my $k (keys %$hashref) {
		$freq{$k} = $hashref->{$k} / $total;
	}
	
	return \%freq;
}

sub entropy_from_counts {
	my ($hashref) = @_;
	
	my $total = 0;
	foreach my $symbol (keys %$hashref) {
		$total += $hashref->{$symbol};
	}
	
	my %freq;
	foreach my $symbol (keys %$hashref) {
		$freq{$symbol} = $hashref->{$symbol} / $total;
	}
	
	return entropy(\%freq);
}

sub entropy {
	my ($hashref) = @_;
	
	my $H = 0;
	foreach my $symbol (keys %$hashref) {
		$H += $hashref->{$symbol} * log($hashref->{$symbol});
	}
	
	return -$H;
}

sub random_dna_gc70 {
	my ($len) = @_;
	
	my $dna = "";
	for (my $i = 0; $i < $len; $i++) {
		my $rnd = rand(1);
		if    ($rnd < 0.15) {$dna .= 'A'}
		elsif ($rnd < 0.30) {$dna .= 'T'}
		elsif ($rnd < 0.65) {$dna .= 'G'}
		else                {$dna .= 'T'}
	}
	return $dna;
}

sub random_dna {
	my ($len) = @_;
	
	my $dna = "";
	for (my $i = 0; $i < $len; $i++) {
		my $rnd = rand(1);
		if    ($rnd < 0.25) {$dna .= 'A'}
		elsif ($rnd < 0.50) {$dna .= 'C'}
		elsif ($rnd < 0.75) {$dna .= 'G'}
		else                {$dna .= 'T'}
	}
	return $dna;
}

sub blank_table {
	my ($mer, $pseudo) = @_;
	
	my $size = 4 ** $mer;
	$pseudo = 0 if not defined $pseudo;
	
	my %H;
	for (my $i = 0; $i < $size; $i++) {
		$H{number_to_dna($i, $mer)} = $pseudo;
	}
	
	return \%H;
}

sub number_to_dna {
	my ($n, $m) = @_;
	
	my $s = "";
	for (my $i = 0; $i < $m; $i++) {
		my $max = 4 ** ($m -$i -1);
		if ($max > $n) {
			$s .= 'A';
		} else {
			my $v = int ($n / $max);
			my $c;
			if    ($v == 0) {$c = 'A'}
			elsif ($v == 1) {$c = 'C'}
			elsif ($v == 2) {$c = 'G'}
			elsif ($v == 3) {$c = 'T'}
			else {die "out of bounds $v"}
			$s .= $c;
			$n -= $v * $max;
		}
	}
	
	return $s;
}

sub dna_to_number {
	my ($seq) = @_;
	
	my $total = 0;
	for (my $i = 1; $i <= length($seq); $i++) {
		my $nt = substr($seq, -$i, 1);
		my $v;
		if    ($nt eq 'A') {$v = 0}
		elsif ($nt eq 'C') {$v = 1}
		elsif ($nt eq 'G') {$v = 2}
		elsif ($nt eq 'T') {$v = 3}
		my $q = $v * (4 ** ($i -1));
		$total += $q;
	}
	
	return $total;
}

sub parse_mformat3 {
	my ($file) = @_;
	
	my %blast;
	my $section = '';
	my $report_num = -1;
	open(IN, $file) or die;
	while (<IN>) {
		if (/^#/) {
			if (/^# (T?BLAST[NXP])/) {
				$blast{program} = $1;
			} elsif (/^# Parameters:/) {
				$section = 'Parameters';
			} elsif (/^# Database: (.+)/) {
				$blast{database} = $1;
				$section = 'Database';
			} elsif (/^#   (\S+)/ and $section eq 'Parameters') {
				my $param = $1;
				if ($param =~ /(\S+)=(\S+)/) {
					$blast{parameter}{$1} = $2;
				} else {
					$blast{parameter}{$param} = 'true';
				}
			}
		} else {
			my ($qid, $sid, $E, $N, $Sprime, $S, $alignlen, $nident,
				$npos, $nmism,  $pcident, $pcpos, $qgaps, $qgaplen, $sgaps,
				$sgaplen, $qframe, $qstart, $qend, $sframe, $sstart, $send,
				$group) = split;
			push @{$blast{result}{$qid}{$sid}},
				{
					expect           => $E,
					hsps             => $N,
					bits             => $Sprime,
					score            => $S,
					align_length     => $alignlen,
					identities       => $nident,
					positives        => $npos,
					percent_identity => $pcident,
					percent_positive => $pcpos,
					query_gaps       => $qgaps,
					query_gaplen     => $qgaplen,
					sbjct_gaps       => $sgaps,
					sbjct_gaplen     => $sgaplen,
					q_frame          => $qframe,
					q_start          => $qstart,
					q_end            => $qend,
					s_frame          => $sframe,
					s_start          => $sstart,
					s_end            => $send,
					group            => $group,
				};
			
		}
	}
	close IN;
	
	return \%blast;	
}


sub cluster_blast {
	my ($blast) = @_;
	
	my %q2s;
	my %s2q;
	foreach my $query (keys %{$blast->{result}}) {
		foreach my $sbjct (keys %{$blast->{result}{$query}}) {
			$q2s{$query}{$sbjct}++;
			$s2q{$sbjct}{$query}++;
		}
	}

	my @sort = sort { keys %{$q2s{$b}} <=> keys % {$q2s{$a}} } keys %q2s;
	my %kill;
	my @cluster;
	
	while (@sort) {
		my $head = shift @sort;
		next if exists $kill{$head};
		
		my %body;
		my @hits = keys %{$q2s{$head}};
		foreach my $hit (@hits) {
			foreach my $match (keys %{$s2q{$hit}}) {
				$kill{$match} = 1;
				$body{$match} = 1;
			}
		}
		
		push @cluster, {
			head => $head,
			body => \%body,
		};
	}

	return \@cluster;
}

sub translate {
	my ($seq) = @_;
	my @aa; # translation kept as a list
	my $translate_length = length($seq) - (length($seq) % 3);
	for(my $i = 0; $ i< $translate_length; $ i+= 3) {
		my $codon = substr($seq, $i, 3);
		if (defined $AA{$codon}) {push @aa, $AA{$codon}}
		else                     {push @aa, 'X'}
	}
	return join("", @aa);
}



1;

###############################################################################

package IK::Span;
use strict; use warnings;

my %FIELD = qw(min 1 max 1 strand 1 type 1 score 1 spans 1 tags 1);

sub new {
	my ($class, %p) = @_;
	my $self = bless {};
	foreach my $key (keys %p) {
		if (exists $FIELD{$key}) {$self->{$key} = $p{$key}}
		else {die "unknown parameter ($key) in IK::Span::new"}
	}
	update($self);
	return $self;
}

sub min    {shift->{min}}
sub max    {shift->{max}}
sub strand {shift->{strand}}
sub type   {shift->{type}}
sub score  {shift->{score}}
sub spans  {shift->{spans}}
sub tags   {shift->{tags}}
sub add_spans {
	my ($self, @span) = @_;
	push @{$self->{spans}}, @span;
	update($self);
}

sub update {
	my $self = shift;
	if (not defined $self->{min})    {$self->{min} = 0}
	if (not defined $self->{max})    {$self->{max} = 0}
	if (not defined $self->{strand}) {$self->{strand} = '.'}
	if (not defined $self->{score})  {$self->{score} = 0}
	if (not defined $self->{spans})  {$self->{spans} = []}
	if (@{$self->{spans}}) {
		my $min = 1e50;
		my $max = 0;
		my $pstrand = 0;
		my $mstrand = 0;
		my $ustrand = 0;
		my $score = 0;
		foreach my $span (@{$self->{spans}}) {
			if ($span->{min} < $min) {$min = $span->{min}}
			if ($span->{max} > $max) {$max = $span->{max}}
			$pstrand++ if $span->{strand} eq '+';
			$mstrand++ if $span->{strand} eq '-';
			$ustrand++ if $span->{strand} eq '.';
			$score += $span->{score};
		}
		$self->{min} = $min;
		$self->{max} = $max;
		if    ($pstrand == @{$self->{spans}}) {$self->{strand} = '+'}
		elsif ($mstrand == @{$self->{spans}}) {$self->{strand} = '-'}
		else                                  {$self->{strand} = '.'}
		$self->{score} = $score / @{$self->{spans}};
	}
}

sub overlap {
	my ($s1, $s2) = @_;
	
	my $coor = 0;
	if ($s1->{min} >= $s2->{min} and $s1->{min} <= $s2->{max}) {$coor = 1}
	if ($s1->{max} >= $s2->{min} and $s1->{max} <= $s2->{max}) {$coor = 1}
	return 0 unless $coor;
	
	if ($s1->{strand} eq $s2->{strand}) {return 1}
	else                                {return -1} # coor but not strand
	return 0;
}

1;
