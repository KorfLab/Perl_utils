package Keith;

# a module of frequently used bioinformatics subroutines used by Keith Bradnam
# Last updated by: $Author$
# Last updated on: $Date$

use strict;


# simple subroutine to just count nucleotide frequencies
# also counts N's separately, and finally works out if
# there is anything other than A,T,C,G or N
sub get_mononucleotides{
        
	# receive sequence, make sure it is all uppercase
    my $seq = uc(shift);
	my $length = length($seq);

   	my ($a,$c,$g,$t,$n,$o) = (0,0,0,0,0,0);
      
    # count mononucleotide frequencies
    $a = ($seq =~ tr/A/A/); 
    $c = ($seq =~ tr/C/C/); 
    $g = ($seq =~ tr/G/G/); 
    $t = ($seq =~ tr/T/T/); 
    $n = ($seq =~ tr/N/N/);

  	# calculate remainder (in case there are other characters present)
	$o = $length - $a - $c - $g -$t - $n;
	
	return($a,$c,$g,$t,$n,$o);
}


# simple subroutine to just count dinucleotide frequencies
# ignores any non-A,C,G,T character
sub get_dinucleotides{
        
	# receive sequence, make sure it is all uppercase
    my $seq = uc(shift);
	my $length = length($seq);
	
	# initialize hash
	my %di;
   	($di{'AA'},$di{'AC'},$di{'AG'},$di{'AT'}) = (0,0,0,0,0);
   	($di{'CA'},$di{'CC'},$di{'CG'},$di{'CT'}) = (0,0,0,0,0);
   	($di{'GA'},$di{'GC'},$di{'GG'},$di{'GT'}) = (0,0,0,0,0);
   	($di{'TA'},$di{'TC'},$di{'TG'},$di{'TT'}) = (0,0,0,0,0);

	# want to keep track of any dinucleotides that have something other than A,C,G,or T.
	my $other=0;

    # count dinucleotide frequencies
 	foreach my $i (0..length($seq)){
    	my $tmp = substr($seq,$i,2);                
        $di{$tmp}++;

		# we might see things like 'NA' or 'CN' so we will keep a count of these as well
		$other++ if ($tmp !~ m/^[ACGT]*$/);
     }
	
	return($di{'AA'},$di{'AC'},$di{'AG'},$di{'AT'},$di{'CA'},$di{'CC'},$di{'CG'},$di{'CT'},$di{'GA'},$di{'GC'},$di{'GG'},$di{'GT'},$di{'TA'},$di{'TC'},$di{'TG'},$di{'TT'},$other);
		
}







1;

__END__

=head1 NAME

Keith;

=head1 DESCRIPTION

Keith is a package that has some basic bioinformatic functions which are mostly applicable to working
with FASTA files

=head1 FUNCTIONS

=head2 get_mononucleotides

Returns a list of counts of A,C,G,T,N, and 'Other' bases in a specified string (presumed to be DNA). 

'Other' is a count of everything which isn't a A,C,G,T or N, and while it potentially
may include other nucleotide ambiguity characters (R,Y,S,W etc.) it can be used
as a crude error checking method (if there are more 'Other' characters than everything else,
it's probably a protein sequence; if there are zero Ts but plenty of 'Other', you may have
an RNA sequence)

Example: 
 # a,c,g, and t are the first 4 characters returned
 my ($a,$c,$g,$t) = get_mononucleotides("aaccacacaggaggattttgaaag"); 

 # but six values are always returned for any string
 my ($a,$c,$g,$t,$n,$other) = get_mononucleotides($sequence); 

=head1 AUTHOR

Keith Bradnam (krbradnam@ucdavis.edu)

=head1 ACKNOWLEDGEMENTS

This software was developed at the Genome Center, UC Davis 
Davis, CA, USA

=head1 COPYRIGHT

Copyright (C) 2008. Keith Brandam. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut
