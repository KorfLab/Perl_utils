#!/usr/bin/perl
#
# kl_shuffler.pl 
#
# a script to calculate the KL distance between two files (A & B) and then:
# 1) calculate the *intra* sequence KL distance for file A (randomly split sequences in file A into two new files)
# 2) calculate the *intra* sequence KL distance for file B
# 3) calculate the distance between a shuffled version of A, and a shuffled version of B
# 4) repeat steps 1-3 N times
# 5) see how many times the intra KL distances exceed the known distance

# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;

die "
kl_shuffler.pl - a program for calculating K-L distance between two FASTA files
and then comparing this distance with K-L distances derived from intra-file comparisons
of sequences in both files.

Usage: <fasta file 1> <fasta file 2> <word size> <number of shuffles>\n" unless @ARGV == 4;


my ($FASTA1,$FASTA2,$WORD,$SHUFFLES) = @ARGV;
die "Please specify a word size of 10 or less\n" if ($WORD >10);



##########################################################################
#
# First of all, read two input files and store sequences in arrays
# These sequences will be used a lot by the rest of the script
#
##########################################################################

my (@master_seqs1,@master_seqs2);

open(FILE,"$FASTA1") || die "Can't open $FASTA1\n";
my $fasta = new FAlite(\*FILE);

# loop through each sequence in target file and add to array
while(my $entry = $fasta->nextEntry) {  
        push(@master_seqs1,uc($entry->seq));
}
close(FILE);


open(FILE,"$FASTA2") || die "Can't open $FASTA2\n";
$fasta = new FAlite(\*FILE);

# loop through each sequence in target file and add to array
while(my $entry = $fasta->nextEntry) {  
        push(@master_seqs2,uc($entry->seq));
}
close(FILE);

# print some warnings if these files contain very few sequences
print "\n";

if(@master_seqs1 < 3){
	print "WARNING: 1st file contains less than 3 sequences, and so shuffling can not be performed\n";
}
elsif(@master_seqs1 <10){
	print "WARNING: 1st file contains fewer than 10 sequences. Please note that shuffling experiements will therefore not be very useful\n";
}
if(@master_seqs2 < 3){
	print "WARNING: 2nd file contains less than 3 sequences, and so shuffling can not be performed\n";
}
elsif(@master_seqs2 <10){
	print "WARNING: 2nd file contains fewer than 10 sequences. Please note that shuffling experiements will therefore not be very useful\n";
}



##########################################################################
#
# Now we can do the first calculation to get the KL distance *BETWEEN* the two files
#
##########################################################################

# need to keep track that this is the first KL distance calculation
my $first_run = "yes";

# will use references for each sequence array
my $comparison_distance = kl(\@master_seqs1,\@master_seqs2);

# can now forget that this was the first run
$first_run = "no";

print "\nInter-sequence K-L distance between file 1 and file 2 using word size $WORD = $comparison_distance\n\n";


# need to count how many times the comparison KL distance is exceeded
my $exceeded = 0;


##########################################################################
#
# Now need to calculate the the KL distance *WITHIN* each of the two files.
#
##########################################################################

# a simple loop to go through sequences from each of the input files

for(my $i=1;$i<3;$i++){
	
	# double check number of sequences in each file - can't shuffle if there is only 1 sequence
	if($i == 1 && @master_seqs1 < 3){
		next;
	}	
	elsif($i == 2 && @master_seqs2 < 3){
		next;
	}
	
	# reset counter
	$exceeded = 0;
	
	# for each run, will want to be able to store randomized arraya of @master_seqs1 & @master_seqs2
	# and to do this we will first make a local copy of these arrays in @seqs
	my (@seqs,@random1,@random2);

	
	# now make lots of pairs of files by randomly selecting sequence from @seqs
	for (my $j=1;$j<=$SHUFFLES;$j++){

		# reset arrays
		@random1 = ();
		@random2 = ();
		
		# make local copy of relevant sequence array
		(@seqs = @master_seqs1) if ($i == 1);
		(@seqs = @master_seqs2) if ($i == 2);
	
		while(@seqs){
			# choose a random position
			my $rand1 = int(rand(1) * @seqs);
						
			# add a sequence to first new random array, and remove sequence from @seqs array 
			push(@random1,$seqs[$rand1]);
			splice(@seqs,$rand1,1);

			# repeat for 2nd randomly chosen sequence (if there is still something in @seqs)
			if(@seqs){
				my $rand2 = int(rand(1) * @seqs);
				push(@random2,$seqs[$rand2]);
				splice(@seqs,$rand2,1);
			}
		}

		# now calculate kl distance from pair of new files
		my $distance = kl(\@random1,\@random2);
		#print "$i $distance\n";

		# has this exceeded first distance? If so, increment counter
		$exceeded++ if ($distance > $comparison_distance);

	}
	print "Inter-sequence distance exceeded $exceeded times out of $SHUFFLES intra-sequence comparisons for file $i\n";
}


##########################################################################
#
# Now need to calculate the KL distance between two files with randomized sequences
#
##########################################################################


# reset counter
$exceeded = 0;

for (my $i=1;$i<=$SHUFFLES;$i++){
	
	# need local copies of arrays for each shuffling run plus two arrays for storing shuffled sequences
	my @seqs1 = @master_seqs1;
	my @seqs2 = @master_seqs2;
	my (@random1,@random2);
	
	# loop through each sequence in input arrays, and randomize sequence before making KL distance calculation

	while(my $seq = shift(@seqs1)) {
        # randomize sequence
        my $random = randomize($seq);
		push(@random1,$random);
		
	}
	while(my $seq = shift(@seqs2)) {
        # randomize sequence
        my $random = randomize($seq);
		push(@random2,$random)
	}

	# now calculate kl distance from pair of new files
	my $distance = kl(\@random1,\@random2);

	# has this exceeded first distance? If so, increment counter
	$exceeded++ if ($distance > $comparison_distance);
}

print "Inter-sequence distance exceeded $exceeded times out of $SHUFFLES when comparing shuffled versions of both files\n\n";


exit(0);

################################################
#
#
#  T H E   S U B R O U T I N E S
#
#
################################################


# the main parent KL distance subroutine
# this gets frequency tables for each array of sequences being compared
# then calculate two KL distances before taking an average distance and converting to bits (nats?)

sub kl{
	my $seq1 = shift;
	my $seq2 = shift;
	
	# get frequencies of words in each sequence file
	my $freq1     = &frequency_table($seq1,$WORD,1);
	my $freq2     = &frequency_table($seq2,$WORD,2);

	# calculate reciprocal K-L distances
	my $distance1 = &kl_distance($freq1,$freq2);
	my $distance2 = &kl_distance($freq2,$freq1);

	# take average of both distances
	my $distance = ($distance1 + $distance2) /2;

	# convert score to bits
	$distance /= log(2);
}




sub frequency_table{
	my ($file,$w,$file_number) = @_;

	# keep track of all words
	my $total_words = 0;

	# keep track of each word, initially set to 1 to avoid zero errors
	my %count = &word_table($w,1);
	
	# loop through each sequence in target file
	foreach my $seq (@{$file}){	
		
		# loop through sequence in windows equal to motif size
	    for (my $i = 0; $i < length($seq)-$w+1; $i++) {
			
			# extract a window of sequence, split it, and place in array    
	    	my $word = substr($seq, $i, $w);

			# check that word is a valid sequence word, count that words and increment total counter
			if (exists $count{$word}){
				$count{$word}++;
				$total_words++;
			}
        }
	}
		
	# check (and warn) if too many words do not exist in the input sequence files
	# it is useful to know if many of the different possible words only exist as pseudocounts
	# Also need to convert counts to frequencies
	my %freq;
	
	# counter for how many words only exist as pseudocounts
	my $only_pseudocounts = 0;
	
	foreach my $word (keys %count){
		$freq{$word} = $count{$word}/$total_words;
		$only_pseudocounts++ if ($count{$word} == 1);
	}

	# now warn if KL distance was based on too many pseudocounts
	# only do this for first kl distance calculation
	if($first_run eq "yes"){
		my $total_keys = keys(%count);
		my $percentage = sprintf("%.0f",$only_pseudocounts/$total_keys *100);
	
		# Only want to warn if a certain proportion of words only exist with counts of 1 (only pseudocounts)
		# The K-L distance may be less meaningful if it based on many comparisons of 1 vs 1 counts
		# The threshold level (5%) is arbitrary but should help to give a clue as to whether using a smaller word size would be more approrpriate

		my $threshold = 5;
	
		if ($only_pseudocounts && ($percentage >= $threshold)){
			print "WARNING: $only_pseudocounts out of $total_keys words ($percentage%) do not exist in file $file_number. Consider using a smaller word size in order to calculate a more reliable K-L distance\n"; 		
		}
	}

	return \%freq;
}

#####################################################
# pseudocount all words
#####################################################

sub word_table{
	my ($w,$pseudocount) = @_;
	my @alphabet = qw(A C G T);
	
	my %table;
	
	# populate tables with pseudocounts, do this for each possible word size
	if($w == 1){
		foreach my $c1 (@alphabet){
			$table{"$c1"} = $pseudocount;
		}
	}
	
	elsif($w == 2){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				$table{"$c1$c2"} = $pseudocount;
			}
		}
	}
	
	elsif($w == 3){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					$table{"$c1$c2$c3"} = $pseudocount;
				}
			}
		}
	}
	
	elsif($w == 4){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						$table{"$c1$c2$c3$c4"} = $pseudocount;
					}
				}
			}
		}
	}
	
	elsif($w == 5){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						foreach my $c5 (@alphabet){
							$table{"$c1$c2$c3$c4$c5"} = $pseudocount;
						}
					}
				}
			}
		}
	}
	
	elsif($w == 6){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						foreach my $c5 (@alphabet){
							foreach my $c6 (@alphabet){
								$table{"$c1$c2$c3$c4$c5$c6"} = $pseudocount;
							}
						}
					}
				}
			}
		}
	}
	
	elsif($w == 7){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						foreach my $c5 (@alphabet){
							foreach my $c6 (@alphabet){
								foreach my $c7 (@alphabet){
									$table{"$c1$c2$c3$c4$c5$c6$c7"} = $pseudocount;
								}
							}
						}
					}
				}
			}
		}
	}
	
	elsif($w == 8){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						foreach my $c5 (@alphabet){
							foreach my $c6 (@alphabet){
								foreach my $c7 (@alphabet){
									foreach my $c8 (@alphabet){
										$table{"$c1$c2$c3$c4$c5$c6$c7$c8"} = $pseudocount;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	elsif($w == 9){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						foreach my $c5 (@alphabet){
							foreach my $c6 (@alphabet){
								foreach my $c7 (@alphabet){
									foreach my $c8 (@alphabet){
										foreach my $c9 (@alphabet){
											$table{"$c1$c2$c3$c4$c5$c6$c7$c8$c9"} = $pseudocount;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	elsif($w == 10){
		foreach my $c1 (@alphabet){
			foreach my $c2 (@alphabet){
				foreach my $c3 (@alphabet){
					foreach my $c4 (@alphabet){
						foreach my $c5 (@alphabet){
							foreach my $c6 (@alphabet){
								foreach my $c7 (@alphabet){
									foreach my $c8 (@alphabet){
										foreach my $c9 (@alphabet){
											foreach my $c10 (@alphabet){
												$table{"$c1$c2$c3$c4$c5$c6$c7$c8$c9$c10"} = $pseudocount;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return %table;	
}

###############################################################


sub kl_distance{
	my ($p,$q) = @_;
	
	my $distance = 0;
	
	foreach my $w (keys %{$p}){
		$distance += $p->{$w} * log ($p->{$w} / $q->{$w});
	}
	return $distance;
}



sub randomize{
    my $seq = shift;
    my $length = length($seq);
    
    # split sequence into an array
    my @seq = split(//,$seq);
    
    # want a random sequence array to store new sequence
    my @random;

    # loop through input sequence, randomly moving k-mers into a new randomized sequence
    while(@seq){
    
		 # choose random coordinate in sequence
		 my $rand = int(rand(1) * @seq);
        
		# add chosen random nt to second array
		push(@random,$seq[$rand]);    
		
	    # and remove from @seq
	    splice(@seq,$rand,1);
    }
    my $random = join('',@random);
    
    return $random;
}

__END__

