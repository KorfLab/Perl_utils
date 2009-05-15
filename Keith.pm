package Keith;

# a module of frequently used bioinformatics subroutines used by Keith Bradnam
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use List::Util qw(sum);

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



# calculate a z test statistic from two sets of numbers.
# will also calculate (and return), the means, standard deviations, 
# and standard errors of the means for both datasets
sub z_test{
	my $a = shift;
	my $b = shift;
	
	# 3rd value can be an optional level of precision
	# defaults to 4 dp
	my $precision = shift;
	$precision = 4 if (!$precision);
	
	# and now all quite simple calculations....
	my $count_a = @$a;
	my $count_b = @$b;
	
	my $sum_a = sum(@$a);
	my $sum_b = sum(@$b);
	
	my $mean_a = sprintf("%.${precision}f",$sum_a/$count_a);
	my $mean_b = sprintf("%.${precision}f",$sum_b/$count_b);
		
	my $std_dev_a = sprintf("%.${precision}f",sqrt(sum(map {($_ - $mean_a) ** 2} @$a) / ($count_a -1)));
	my $std_dev_b = sprintf("%.${precision}f",sqrt(sum(map {($_ - $mean_b) ** 2} @$b) / ($count_b -1)));

    my $std_err_a = sprintf("%.${precision}f",$std_dev_a/sqrt($count_a));
    my $std_err_b = sprintf("%.${precision}f",$std_dev_b/sqrt($count_b));

	# ... leading up to the final z statistic
    my $z = sprintf("%.${precision}f",($mean_a - $mean_b)/sqrt(($std_dev_a**2/$count_a)+($std_dev_b**2/$count_b)));

	return($z,$mean_a,$mean_b,$std_dev_a,$std_dev_b,$std_err_a,$std_err_b);
} 

#adds a new line character every 60 bases  

sub tidy_seq{
    my ($seq) = @_;
    $seq =~ s/[\s\n]//g;
    $seq =~ tr/a-z/A-Z/;
    
    my ($output_seq) = "";
    my (@seq2) = "";
    my ($start,$end);

    @seq2 = split(//,$seq);
    my $length = @seq2;
    my ($to_add) = int($length/60);

    $end = $start= 0;

    foreach (1..$to_add){
        $output_seq .= substr($seq,$start,60);
        $output_seq .= "\n";
        $start += 60;
    }
    $output_seq .= substr($seq,$start);
    return ($output_seq);
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




=head2 get_dinucleotides

Returns a list of counts of all 16 dinucleotides, and 'Other' for any specified string (presumed to be DNA). 

'Other' is a count of any dinucleotide which includes something other than A,C,T, or G. The percentage
calculations exclude these 'Other' dinucleotides. 

Example: 
 my (@dinucs) = get_dinucleotides("aaccacacaggaggattttgaaag"); 


=head2 z_test

If you have two lists of numbers (in separate arrays), send them to the z_test function
and it will calculate a z test statistic which will let you know how different the means
of the two samples are (and whether the difference is significant).

It will also (if you want) return the information used to calculate the z statistic (means,
standard deviations, and standard errors)

The function defaults to a level of 4 decimal places for precision, but you change this with 
the third argument to the function

Example (increasing level of precision to 5 dp): 
my ($z,$mean_a,$mean_b,$std_dev_a,$std_dev_b,$std_err_a,$std_err_b) = Keith::z_test(\@a,\@b,5);


=head2 tidy_seq

Simply pass subroutine a DNA/RNA/protein sequence (just the sequence, not the header) as a 
string and subroutine will tidy up sequence to be all upper case and just 60 characters per line.



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
