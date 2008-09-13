#!/usr/bin/perl
#
# kl_randomizer.pl 
#
# a script to calculate the KL distance between two files (A & B) and then:
# 1) calculate the distance between a randomized version of A, and a randomized version of B
# 2) repeat step 1 N times
# 3) see how many times the intra KL distances exceed the known distance
# 4) quantify the distribution of KL distances from the randomized sequences, and assign a probability to the real distance

# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;
use List::Util qw(sum);

die "
kl_randomizer.pl - a program for calculating K-L distance between two FASTA files
and then comparing this distance with K-L distances derived from randomized
versions of both files (preseving mononucleotide frequencies).

Usage: <fasta file 1> <fasta file 2> <word size> <number of randomizations>\n" unless @ARGV == 4;


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





##########################################################################
#
# Now we can do the first calculation to get the KL distance *BETWEEN* the two files
#
##########################################################################

# need to keep track that this is the first KL distance calculation
my $first_run = "yes";

# will use references for each sequence array
my $comparison_distance = sprintf("%.6f",kl(\@master_seqs1,\@master_seqs2));

# can now forget that this was the first run
$first_run = "no";

print "\nInter-sequence K-L distance between file 1 and file 2 using word size $WORD = $comparison_distance\n\n";




##########################################################################
#
# Now need to calculate the KL distance between two files with randomized sequences
#
##########################################################################

print "Randomzing sequences $SHUFFLES times:\n";

# need to count how many times the comparison KL distance is exceeded
my $exceeded = 0;

# store all distances 
my @distances;

# turn on flushing of output buffer
$| = 1;

for (my $i=1;$i<=$SHUFFLES;$i++){
	if($i == 1){
		print "1";
	}
	elsif($i % 10 ==0){
		print "$i";
	}
	elsif($i == $SHUFFLES){
		print "$SHUFFLES\n";
	}
	else{
		print ".";
	}
	
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

	# now calculate kl distance from pair of new files, and add to array
	my $distance = kl(\@random1,\@random2);
	push(@distances,$distance);

	# has this exceeded first distance? If so, increment counter
	$exceeded++ if ($distance > $comparison_distance);
}
print "\n\n";

print "Inter-sequence distance exceeded $exceeded times out of $SHUFFLES when comparing randomized versions of both files\n\n";


# now calculate a mean and standard deviation of all the individual KL shuffled distances
my $n = @distances;
my $mean = sprintf("%.6f",sum(@distances)/$n);
my $std_dev = sprintf("%.6f",sqrt(sum(map {($_ - $mean) ** 2} @distances) / ($n-1)));
my $z = sprintf("%.2f",($comparison_distance-$mean)/$std_dev);

# display some probability info
my $p = "N.S.";
$p = "<0.05"    if ($z >= 1.96);
$p = "<0.01"    if ($z >= 2.58);
$p = "<0.001"   if ($z >= 3.29);
$p = "<0.0001"  if ($z >= 3.89);
$p = "<0.00001" if ($z >= 4.41);

print "N = $n\n";
print "Mean KL distance = $mean\n";
print "Standard deviation = $std_dev\n";
print "Z score = $z\n";
print "Probability = $p\n";

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

