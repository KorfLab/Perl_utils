#!/usr/bin/perl -w
#
# shuffle_seq.pl
#
# A script to take a file of real worm introns and then write a new file of random sequences
# with the same length as the real introns, but with nucleotide frequencies based on real
# intron nucleotide frequencies
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;


########################
# Command line options
########################

my ($file,$word,$output) = @ARGV;

# sanity checks
die "Usage: shuffle_seq.pl <fasta file> <word size to shuffle> <number of random sequences to create>\n" if (@ARGV != 3);
die "Enter a positive integer for word size and number of output files to create\n" if ($word < 1 || $output < 1);

print "\nCreating $output random sequences from $file, using $word letters as the length unit for shuffling...\n";

open(IN,"<$file") || die "Couldn't open $file\n";

my $fasta = new FAlite(\*IN);

# loop through each sequence in input file

while(my $entry = $fasta->nextEntry) {
    my $header = $entry->def;

	# grab sequence
    my $seq = uc($entry->seq);


	# loop through each possible random output file that we will create
	for(my $i = 1;$i<=$output;$i++){
		open(OUT,">>$file.r$i") || die "Couldn't write output file\n";

	    print OUT "$header\n";

	    # randomize sequence
		my $random = randomize($seq);

		# add newline every 60 characters
		print OUT "$random\n";

		close(OUT) || die "Couldn't close output file\n";

	}

}

close(IN) || die "Couldn't close $file\n";

exit(0);



#########################################################################

sub randomize{
	my $seq = shift;
	my $length = length($seq);
	
	# if word size is bigger or equal to sequence size, then we are in trouble
	if($word >= $length){
		die "Your word size for shuffling ($word nt), is as long (or longer) than one of the sequences in your input file. Please choose a smaller value.\n\n";
	}

	# split sequence into an array
	my @seq = split(//,$seq);

  	
	# want a random sequence array to store new sequence
	my @random;

	# loop through input sequence, randomly moving k-mers into a new randomized sequence
	while(@seq){
		
	
		# what is the maximum position that we could select as a start coordinate to extract sequence
		# this depends on word size (e.g. a word size of 9 from a sequence of 10 nt only has a start coordinate of 0 or 1)
		my $max_pos = @seq - $word + 1;
		
		# choose random coordinate in sequence
		my $rand = int(rand(1) * $max_pos);
			
        # add selected random oligo to second array
		# if the remaining sequence is less than the word size, then we simply just add everything that is left in @seq
		if (@seq < $word){
			 push(@random,@seq);    
			# and remove from @seq
			splice(@seq,$rand,$word);
			
		}
		# otherwise, only work with current k-mer
		else{
			push(@random,@seq[$rand..$rand+$word-1]);  
			# and remove from @seq
			splice(@seq,$rand,$word);
		}
	}
	my $random = join('',@random);
	
	# tidy sequence before returning it (i.e. add newline every 60 characters)
	$random = tidy_seq($random);

	return $random;
}




sub tidy_seq{
    # adds a new line character every 60 bases  
    my $seq = shift;
    $seq =~ s/[\s\n]//g;
    $seq =~ tr/a-z/A-Z/;

    # how many newlines to add?
    my $to_add = int(length($seq)/60);
    
    my $output_seq;     
    my ($start,$end);
    $end = $start= 0;

    foreach (1..$to_add){
        $output_seq .= substr($seq,$start,60);
        $output_seq .= "\n";
        $start += 60;
    }
    # add remaining sequence (if it exists)
    $output_seq .= substr($seq,$start);

    return ($output_seq);
}



__END__