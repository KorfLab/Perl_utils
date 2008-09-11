#!/usr/bin/perl
#
# kl distance calculator.pl 
#
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict;
use warnings;
use FAlite;

die "
kl_distance.pl - a program for calculating K-L distance between two FASTA files
Usage: <fasta file 1> <fasta file 2> <word size>\n" unless @ARGV == 3;


my ($FASTA1,$FASTA2,$WORD) = @ARGV;

die "Please specify a word size of 10 or less\n" if ($WORD >10);


# get frequencies of words in each sequence file
my $freq1     = &frequency_table($FASTA1,$WORD);
my $freq2     = &frequency_table($FASTA2,$WORD);


# calculate reciprocal K-L distances
my $distance1 = &kl_distance($freq1,$freq2);
my $distance2 = &kl_distance($freq2,$freq1);

# take average of both distances
my $distance = ($distance1 + $distance2) /2;

# convert score to bits
$distance /= log(2);

print "\nK-L distance = $distance\n\n";



exit(0);

################################################
#
#
#   T H E   S U B R O U T I N E S
#
#
################################################

sub frequency_table{
	my ($file,$w) = @_;

	# keep track of all words
	my $total_words = 0;

	# open file to read sequence and break into words
	# then count work frequencies to return values

	open(FILE,"$file") || die "Can't open $file\n";

	my $fasta = new FAlite(\*FILE);

	# keep track of each word, initially set to 1 to avoid zero errors
	my %count = &word_table($w,1);
	
	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry) {  
	        my $seq = uc($entry->seq);
			
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
	close(FILE);
	

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

	my $total_keys = keys(%count);
	my $percentage = sprintf("%.0f",$only_pseudocounts/$total_keys *100);
	
	# Only want to warn if a certain proportion of words only exist with counts of 1 (only pseudocounts)
	# The K-L distance may be less meaningful if it based on many comparisons of 1 vs 1 counts
	# The threshold level (5%) is arbitrary but should help to give a clue as to whether using a smaller word size would be more approrpriate

	my $threshold = 5;
	
	if ($only_pseudocounts && ($percentage >= $threshold)){
		print "\nWARNING: in file $file, $only_pseudocounts out of a possible $total_keys words ($percentage%) do not exist at all\n";
		print "Maybe consider using a smaller word size in order to calculate a more reliable K-L distance\n"; 		
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



__END__

