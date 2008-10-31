use strict;
use warnings;
use Storable;
use FAlite;
#FIX: SCORE SEQUENCE FUNCTION
#FIX:  HTML REPORT
#TODO:  IMPLEMENT SVG
#TODO:  IMPLEMENT R
package Motif;

################################################################################
#
#   new
#
#   Description:  Create a new Motif class
#
# 
#

sub new{
    my $class = shift;
    my $self = bless {}, $class;
    return $self;
}

###############################################################################
#
#   _generate_Consensus:
#
#   Description:     generates a consensus sequene from Position Probability Matrix
#
#   Work: get function working and use to fill in consensus
#
sub _generate_Consensus{
    my ($self,$alphabet)=@_;
    my @PPM=@{$self->{POSITION_PROBABILITY_MATRIX}};
    my @Consensus;
##                  4: [ACGT]
##                  8: [ACGTacgt]
##                  11: [ACGTRYMKWSN]
##                  15: [ACGTRYMKWSNBDHV]
##                  21: [ACGTRYMKWSNacgtrymkws]
##                  25: [ACGTRYMKWSBDHVNacgtrymkws]
##                  29: [ACGTRYMKWSBDHVNacgtrymkwsbdhv]

#Hash containing probabilities for the particular ambiguous nucleotide
    my %Alpha4 = (	A	=>	[1,0,0,0], ##[A,C,G,T]
                	C	=>	[0,1,0,0],
                	G	=>	[0,0,1,0],
                	T	=>	[0,0,0,1]);
    my %Alpha11 = (	R	=>	[0,0.5,0.5,0],
        		Y	=>	[0.5,0,0,0.5],
        		M	=>	[0.5,0.5,0,0],
        		K	=>	[0,0,0.5,0.5],
        		W	=>	[0.5,0,0,0.5],
        		S	=>	[0,0.5,0.5,0],
        		N	=>	[.25,.25,.25,.25]);
    my %Alpha15 = (	B	=>	[0,(1/3),(1/3),(1/3)],
        		D	=>	[(1/3),0,(1/3),(1/3)],
        		H	=>	[(1/3),(1/3),0,(1/3)],
        		V	=>	[(1/3),(1/3),(1/3),0]);
    my %Alpha8 = (	a	=>	[0.75,(1/12),(1/12),(1/12)],
        		c	=>	[(1/12),0.75,(1/12),(1/12)],
        		g	=>	[(1/12),(1/12),0.75,(1/12)],
        		t	=>	[(1/12),(1/12),(1/12),0.75]),
    my %Alpha21= (	r	=>	[0.4175,0.0825,0.4175,0.0825],
        		y	=>	[0.0825,0.4175,0.0825,0.4175],
        		m	=>	[0.4175,0.4175,0.0825,0.0825],
        		k	=>	[0.0825,0.0825,0.4175,0.4175],
        		w	=>	[0.4175,0.0825,0.0825,0.4175],
        		s	=>	[0.0825,0.4175,0.4175,0.0825]);
    my %Alpha29 = (	b	=>	[0.125,(7/24),(7/24),(7/24)],
        		d	=>	[(7/24),0.125,(7/24),(7/24)],
        		h	=>	[(7/24),(7/24),0.125,(7/24)],
        		v	=>	[(7/24),(7/24),(7/24),0.125]);
    my %alphabet;

##Append Alphabet hashes to make appropriate hash with corresponding letters
##Appropiate alphabet sizes are: 4, 8, 11, 15, 21, 25, 29
##                  4: [ACGT]
##                  8: [ACGTacgt]
##                  11: [ACGTRYMKWSN]
##                  15: [ACGTRYMKWSNBDHV]
##                  21: [ACGTRYMKWSNacgtrymkws]
##                  25: [ACGTRYMKWSBDHVNacgtrymkws]
##                  29: [ACGTRYMKWSBDHVNacgtrymkwsbdhv]
    if ($alphabet==4)		{%alphabet= %Alpha4;}
    elsif ($alphabet==8)	{%alphabet= (%Alpha4,%Alpha8)}
    elsif ($alphabet==11)	{%alphabet= (%Alpha4,%Alpha11)}
    elsif ($alphabet==15)	{%alphabet= (%Alpha4,%Alpha11,%Alpha15)}
    elsif ($alphabet==21)	{%alphabet= (%Alpha4,%Alpha8,%Alpha11,%Alpha21)}
    elsif ($alphabet==25)	{%alphabet= (%Alpha4,%Alpha8,%Alpha11,%Alpha15,%Alpha21)}
    elsif ($alphabet==29)	{%alphabet= (%Alpha4,%Alpha8,%Alpha11,%Alpha15,%Alpha21,%Alpha29)}
    else {die "Non-supported word size";}

## Converts Counts stored in @Count to appropriate alphabet by measuring K-L distance between Distributions
    foreach my $position (@PPM){  #For each Array in @PPM
            my @KLDistance;  #best K-L distance so far
	    my @Sample_Distr=@{$position};  #removed psuedocount command.....Need to implement it.
            foreach my $keys (keys %alphabet){  #For each of the letters in the alphabet mesure K-L distance
                    #my @Actual = [@{$key}];              #Put first array from the Samples hash into array Actual
                    my @Consensus_distr=@{$alphabet{$keys}};
		    #@values=[pseudocount(@values)];
		    #my $Measure_distance = KL(@Actual,@values,$count);  #Measure distance passes KL subroutine the sample array and F
                    my $KL_Dist;
                    for (my $i=0;$i<4;$i++){
                        if ($Sample_Distr[$i]==0){
                            $Sample_Distr[$i]=0.000000000001;
                        }
                        if ($Consensus_distr[$i]==0){
                            $Consensus_distr[$i]=0.000000000001;
                        }
                        $KL_Dist += $Sample_Distr[$i]*(log($Sample_Distr[$i]/$Consensus_distr[$i])/log(2))+$Consensus_distr[$i]*(log($Consensus_distr[$i]/$Sample_Distr[$i])/log(2));
                    }
		    if (!defined $KLDistance[0]){       #If the K-L Distance is undefined then save it in @KLDistance or else if it's less than what is stored in KLDistance replace it.
                            @KLDistance = ($keys,$KL_Dist)}
                    elsif ($KL_Dist<$KLDistance[1]){
                            @KLDistance = ($keys,$KL_Dist)}
                    else {next;}
            }
            push @Consensus, [@KLDistance];  #The letter and distance that is the smallest distance between distributions are pushed into String
    }
    return @Consensus;
}

#----------------------------------------------------------------------------------------------------------------------------------------------------
package Motifset;
require Storable;
require Cwd;

################################################################################
#
#   new()
#
#   Description:  creates new Motifset class, with option to import Motifs
#                 automatically.
#
#   Usage:  my $new_motifset=Motifset->new();
#           my $new_motifset=Motifset->new("TRANSFAC"=>"Transfac.txt","WEEDER"=>"filename.wee");  import single weeder and transfac file
#           my $new_motifset=Motifset->new("TRANSFAC"=>["Transfac.txt","Transfac.txt"], "WEEDER"=>"filename.wee"); import multiple tranfac matrices and weeder file
#           my $new_motifset=Motifset->new("JASPAR"=>{"ANNOTATION"=>"MATRIX_ANNOTATION.txt","DATA"=>"MATRIX_DATA.txt"}); import Jaspar database
#           my $new_motifset=Motifset->new("JASPAR"=>[{"ANNOTATION"=>"MATRIX_ANNOTATION.txt","DATA"=>"MATRIX_DATA.txt"},{"ANNOTATION"=>"MATRIX_ANNOTATION.txt","DATA"=>"MATRIX_DATA.txt"}]);  #import multiple Jaspar database files
#
#
sub new{
    my %default=("TRANSFAC"=>undef,
                 "WEEDER"=>undef,
                 "NESTEDMICA"=>undef,
                 "MEME"=>undef,
                 "JASPAR"=>undef);
    my $class = shift;
    my %arg=(%default,@_);
    my $self = bless {}, $class;

    if (!defined @_){
        return $self;}
    else {
        foreach my $motif_finder (keys %arg){
            if (!defined $arg{$motif_finder}){
                next;
            }
            else {
                if (ref $arg{$motif_finder} eq "ARRAY") {
                    foreach my $file (@{$arg{$motif_finder}}){
                        $self=  ($motif_finder eq "TRANSFAC") ? transfac_import($self,$file):
                                ($motif_finder eq "WEEDER") ? weeder_import($self,$file):
                                ($motif_finder eq "NESTEDMICA") ? nestedmica_import($self,$file):
                                ($motif_finder eq "MEME") ? meme_import($self,$file):
                                ($motif_finder eq "JASPAR") ? jaspar_import($self,%$file):
                                die "Unknown import parameters $motif_finder and $file";
                    }
                }
                else {
                    $self=  ($motif_finder eq "TRANSFAC") ? transfac_import($self,$arg{$motif_finder}):
                            ($motif_finder eq "WEEDER") ? weeder_import($self,$arg{$motif_finder}):
                            ($motif_finder eq "NESTEDMICA") ? nestedmica_import($self,$arg{$motif_finder}):
                            ($motif_finder eq "MEME") ? meme_import($self,$arg{$motif_finder}):
                            ($motif_finder eq "JASPAR") ? jaspar_import($self,%{$arg{$motif_finder}}):
                            die "Unknown import parameters $motif_finder and $arg{$motif_finder}";
                }
            }
        }
    }
    return $self;
}


################################################################################
#
#   query()
#
#   Description: Given a string or hash of parameters, searches and returns
#               motifset with motifs that match.   If no mofifs then
#
#
#
#
sub query{
    my %categories=("ACCESSION"=>undef,
                    "ANY"=>undef,
                    "ALPHABET"=>undef,
                    "BEST_OCCURRENCES_SCORING_MATRIX"=>undef,
                    "BACKGROUND_FREQUENCY"=>undef,
                    "BINDING_SITES_UNDERLYING_MATRIX"=>undef,
                    "COMMAND_LINE"=>undef,
                    "COMMENTS"=>undef,
                    "CONSENSUS_SEQUENCE"=>undef,
                    "COPYRIGHT"=>undef,
                    "DATAFILE"=>undef,
                    "DATE_AUTHOR_UPDATE_DATES"=>undef,
                    "FILE"=>undef,
                    "FILENAME"=>undef,
                    "IDENTIFIER"=>undef,
                    "JOB"=>undef,
                    "LIST_OF_LINKED_FACTOR_ENTRIES"=>undef,
                    "MEDLINE_ID"=>undef,
                    "MOTIF_JOB"=>undef,
                    "MOTIF_LOCATIONS"=>undef,
                    "MOTIF_METHOD"=>undef,
                    "MOTIF_SEQUENCES"=>undef,
                    "NAME_OF_BINDING_SITE"=>undef,
                    "ORGANISM"=>undef,
                    "POSITION_PROBABILITY_MATRIX"=>undef,
                    "POSITION_SCORING_MATRIX"=>undef,
                    "REFERENCE_#"=>undef,
                    "REFERENCE_AUTHORS"=>undef,
                    "REFERENCE_DATA"=>undef,
                    "REFERENCE_TITLE"=>undef,
                    "REGULAR_EXPRESSION"=>undef,
                    "SEQUENCES"=>undef,
                    "SHORT_FACTOR_DESCRIPTION"=>undef,
                    "STATISTICAL_BASIS"=>undef,
                    "THRESHOLD"=>undef,
                    "VERSION"=>undef);

    my ($self,$motifset,$searchterm,$category)=@_;
    my %motif_accessions;

    if (!defined $category){
        $category="ANY"
    }
    
    foreach my $motif (@{$self->{MOTIFS}}){
        $motif_accessions{$motif->{ACCESSION}}=1;
    }

    #check for compatible categories
    if ($category ne "ANY"){
        if (!exists $categories{$category}){
            $category="ANY"
        }
    }

    #if any then search whole motif until true
    if ($category eq "ANY"){
        foreach my $motif (@{$motifset->{MOTIFS}}) {
            foreach my $key (keys %$motif){
                if (main::_query($motif->{$key},$searchterm)){
                    if (!exists $motif_accessions{$motif->{ACCESSION}}){
                        push @{$self->{MOTIFS}},$motif;
                    }
                    last;
                }

            }
        }
    }

    #else search only the specified category
    if ($category ne "ANY"){
        foreach my $motif (@{$motifset->{MOTIFS}}) {
            if (exists $motif->{$category}){
                if (main::_query($motif->{$category}, $searchterm)){
                    push @{$self->{MOTIFS}},$motif;
                }
            }
            else {
                next;
            }

        }
    }
    return $self;
}




############################ Not finished
#
#   print_summary
#
#   Description:  Prints Summary of Motifs
#   Non-verbose (default)
#
#   Flags:  Default, -v (verbose), -s (accession numbers only)
#   verbose -v  prints everything in Motifset
#   Default:  Prints (Accession, File used to import data, Position Scoring Matrix, Position Probability Matrix, Sequences used to discover Motif)
#
#   Work:  Need to standardize output of the function and improve flag handling
#          Remove accessions
#
#

sub print_summary{
    my ($self,$flag)=@_;
    $|=1;
    foreach my $key (sort keys %$self){
            if (ref $self->{$key} eq "Motif"){
                print "-----------------------------------------------\nMotif: $key\n";
                if ( uc($flag) eq "-S"){
                    next;}
                print "-----------------------------------------------\n";
                print "ACCESSION: " . $key . "\n";
                if (exists $self->{$key}->{IDENTIFIER}){
                    print "IDENTIFIER: " . $self->{$key}->{IDENTIFIER}->[0] . "\n";
                }
                print "Import File: " . $self->{$key}->{FILE}->[0] . "\n";
                if (exists $self->{$key}->{POSITION_SCORING_MATRIX}){
                    print "\n-------------------------\nPosition Scoring Matrix:\n-------------------------\n";
                    printf "%-5s%-5s%-5s%-5s\n","A","C","G","T";
                    foreach my $position (@{$self->{$key}->{POSITION_SCORING_MATRIX}}){
                        printf "%-5d%-5d%-5d%-5d\n",@$position;
                    }
                }
                if (exists $self->{$key}->{POSITION_PROBABILITY_MATRIX}){
                    print "\n-------------------------\nPosition Probability Matrix:\n-------------------------\n";
                    printf "%-5s%-5s%-5s%-5s\n","A","C","G","T";
                    foreach my $position (@{$self->{$key}->{POSITION_PROBABILITY_MATRIX}}){
                        printf "%-5.2f%-5.2f%-5.2f%-5.2f \n",@$position;
                    }
                }
                if(uc($flag) eq "-V"){
                    my @prompts=("IDENTIFIER",
                             "Motif_JOB",
                             "DATE_AUTHOR_UPDATE_DATES",
                             "COPYRIGHT","NAME_OF_BINDING_SITE",
                             "SHORT_FACTOR_DESCRIPTION",
                             "LIST_OF_LINKED_FACTOR_ENTRIES",
                             "STATISTICAL_BASIS",
                             "BINDING_SITES_UNDERLYING_MATRIX",
                             "COMMENTS",
                             "MEDLINE_ID",
                             "REFERENCE_#",
                             "REFERENCE_AUTHORS",
                             "REFERENCE_TITLE",
                             "REFERENCE_DATA",
                             "Motif_SEQUENCES",
                             "SEQUENCES",
                             "Motif_METHOD",
                             "COMMAND_LINE",
                             "DATAFILE",
                             "ALPHABET",
                             "REGULAR_EXPRESSION",
                             "CONSENSUS_SEQUENCE",
                             "VERSION",
                             "THRESHOLD",
                             "ORGANISM");
                    foreach my$prompt (@prompts){
                        if (exists $self->{$key}->{$prompt}){
                            print "\n-------------------------\n$prompt\n-------------------------\n" ;
                            print join "\n", @{$self->{$key}->{$prompt}};
                        }
                        if (exists $self->{$key}->{Motif_LOCATIONS}){
                            ##fill in location print
                        }
                    }
                print "\n\n";
            }
        }
    }
    print "\nEnd of Print Summary\n\n\n\n\n"
}

##############################################################################
#
#   delete
#
#
#   Description:  delete_Motif: deletes a Motif from the Motifset
#  
sub delete{
    my ($self,$delete)=@_;

    if (!defined $delete) {
        my @list;
        my $input;
        foreach my $motif (@{$self->{MOTIFS}}){
            push @list, $motif->{ACCESSION};
        }

        LIST:
        my $n=0;
        foreach (@list){
            print "$n\) $_\n";
            $n++;
        }
        print "Enter the Motif that you would like to delete( Example: 3):";
        $input=<STDIN>;
        if (!defined $input) {
            goto LIST;
        }
        chomp $input;
        if (exists $self->{MOTIFS}->[$input]){
            splice(@{$self->{MOTIFS}},$input,1);
            print "$list[$input] deleted from Motifset\n";
        }
        else {
            print [$input] . " doesn't exist as an option.";
        }

    }
    elsif (exists $self->{MOTIFS}->[$delete]) {
        splice(@{$self->{MOTIFS}},$delete,1);
    }
    elsif (!exists $self->{MOTIFS}->[$delete]){
        print "$delete doesn't exist.";
    }
    else {
        print "Motif doesn't contain $delete\n";
    }

    return $self;
}


##############################################################################
#
#   combine
#
#   Description: takes 2 Motifsets and combines them into one Motifset
#  
sub combine{
    my ($self1,$self2)=@_;

    if (exists $self1->{Distance_Score}){
        delete $self1->{Distance_Score};
    }

    my @motif2=@{$self2->{MOTIFS}};
    push @{$self1->{MOTIFS}},@motif2;

    return $self1;
}


###############################################################################
#
#   import()
#
#   Description: imports a Motif into a Motifset
#
#   Work: Remove dependence upon accession

sub import{
    my ($self,$Motif)=@_;
    if (ref $Motif eq "Motifset"){
        $Motif=main::_query_user($Motif);
    }
    my $accession=$Motif->{ACCESSION};
    my $header=substr($accession,0,2);
    if (exists $self->{$accession}){
        print "The accession $accession conflicts with another accession.\n";
        $self->_import_Motif($Motif);
        print "Accession changed from $accession to $Motif->{ACCESSION}\n";
    }
    else {
        $self->{$accession}=$Motif;
    }
    return $self;
}



###############################################################################
#
#   subset
#
#   Description:    Create subset of Motifset; Inputs: reference to array or
#                   user queried input.
#
#   
sub subset{
    my ($self,$motifset,$subset)=@_;
    if (!defined $motifset) {
        die "No Motifset handed to function to take subset of\n Correct usage\:  \$Motifset->subset(\$meme) \n" ;
    }
    my @subset;
    my @list;
    foreach my $motif (@{$self->{MOTIFS}}){
            push @list, $motif->{ACCESSION};
        }

    if (defined $subset) {
        @subset=@$subset;
        @subset=sort { $a <=> $b } @subset;}

    else{
        my $n=0;
        foreach (@list){
            print "$n\) $_\n";
            $n++;
        }
        print "Enter the Motifs that you would like to create subset( Example: 1,2,3,4):";
        my $input=<STDIN>;
        chomp $input;
        @subset=split /,/,$input;
        @subset=sort { $a<=>$b } @subset;
    }

    for (my $i=1;$i<=scalar @list;$i++){
        if (!defined $subset[0]){
            last;}
        if ($i==$subset[0]){
            shift @subset;
            push @{$self->{MOTIFS}}, $motifset->{MOTIFS}->[$i-1];
            next;}
        else { next;}
    }
    return $self;
}

###############################################################################
#
#   save()
#
#   Description:    save object to storable format under filename
#
#   Work:  change name and function to reflect import;  Possibly depreciate func
#          and incorate into transfac-like output file
sub save{
    my ($self,$filename)=@_;
    Storable::store ($self, "$filename.mts");
    print "Motifset saved to $filename\n";
}

###############################################################################
#
#   retrieve
#
#   Description:  import motifset from storable filename
#
#   Work:  change name to reflect import function; Possibly depreciate and
#           make function for import/export to transfac-like file

sub retrieve{
    my ($class,$filename)=@_;
    my $self=Storable::retrieve("$filename.mts");
    print "Motifset retrieved from $filename\n";
    return $self;
}


#############################################################################
#
#   score2motifsets()
#
#   Description:  Compare motifsets by calculating distance
#
#

sub score2motifsets{
    my %defaults=("database"=>undef,
                  "score"=>"NM",
                  "align"=>"global",
                  "background"=>{'A'=>0.25, 'C'=>0.25, 'G'=>0.25, 'T'=>0.25},
                  "counts"=>undef,
                  "trim"=>4,
                  "output_file"=>undef);
    my $query=shift;


    my %arg=(%defaults,@_);

    my $database=$arg{database};
    my $scoring_method=$arg{score};
    my $alignment_method=$arg{align};
    my $trim=$arg{trim};

    $arg{background}=($defaults{background} == $arg{background}) ? $arg{background}:  #if seq_freq then get single nucleotide frequency table
                    $arg{background}->{1};

    my @background;
    foreach my $letter (sort keys %{$arg{background}}){
        push @background, ${$arg{background}}{$letter};
    }

    my $align_method=   ($alignment_method eq "global") ? "global" :
                        ($alignment_method eq "gapped") ? "gapped"   :
                        ($alignment_method eq "trimmed") ? "trimmed"  :
                        die "Unsupported scoring method";

    my $score_method=   ($scoring_method eq "NM") ? "NestedMica" :
                        ($scoring_method eq "KLD") ? "Kullback-Leibler Distance" :
                        ($scoring_method eq "PCC") ? "Pearson's Correlation Coefficient" :
                        ($scoring_method eq "SW") ? "Sandelin/Wasserman" :
                        ($scoring_method eq "ALLR") ? "Average Log-likelihood Ratio" :
                        ($scoring_method eq "ALLR2") ? "Average Log-likelihood Ratio with limit of 2 imposed":
                        ($scoring_method eq "PCT") ?  "Pearson's Chi-Square Test" :
                        ($scoring_method eq "FIET") ? "Fisher-Irwin Exact Test" :
                        ($scoring_method eq "ED") ? "Euclidean Distance" :
                        die "Unsupported scoring method";

    my @score_table;
    for (my $i=0;$i<scalar(@{$query->{MOTIFS}});$i++){
        my $query_label=$query->{MOTIFS}->[$i]->{ACCESSION};
        my $database_label;
        my $score=0;
        my @score_line;
        my $j=$i+1;
        my @p_ppm1=@{$query->{MOTIFS}->[$i]->{POSITION_PROBABILITY_MATRIX}};
        my @counts;

        #Find pseudocounted counts from scoring matrix, because these functions require counts.
        #First checks if the function requiring counts are called
        #Second checks to see if POSITION_SCORING_MATRIX is available for function.
        #If the PSM isn't available and counts weren't handed to function then die
        #If counts were handed to function then it uses these counts.
        if ((($scoring_method eq "ALLR") || ($scoring_method eq "ALLR2")||
             ($scoring_method eq "FIET") || ($scoring_method eq "PCT"))){
            if (exists $query->{MOTIFS}->[$i]->{POSITION_SCORING_MATRIX}) {
                my $sum+=$_ foreach @{$query->{MOTIFS}->[$i]->{POSITION_SCORING_MATRIX}->[0]};
                $sum = $sum + ($sum)**0.5;
                push @counts,$sum;
            }
            elsif ((!exists $query->{MOTIFS}->[$i]->{POSITION_SCORING_MATRIX}) && ($arg{counts}==$defaults{counts})) {
                die "Counts missing for motifs created in NestedMica.  $score_method requires counts.\n   Please pass parameters:  counts => x, where x is number of counts\n ";
            }
            else {
                push @counts, $arg{counts};
            }
        }

        for (my $j=0;$j<scalar(@{$database->{MOTIFS}});$j++){
            my @p_ppm2=@{$database->{MOTIFS}->[$j]->{POSITION_PROBABILITY_MATRIX}};
            $database_label=$database->{MOTIFS}->[$j]->{ACCESSION};


            #Find pseudocounted counts from scoring matrix, because these functions require counts.
            #First checks if the function requiring counts are called
            #Second checks to see if POSITION_SCORING_MATRIX is available for function.
            #If the PSM isn't available and counts weren't handed to function then die
            #If counts were handed to function then it uses these counts.
            if ((($scoring_method eq "ALLR") || ($scoring_method eq "ALLR2")||
                 ($scoring_method eq "FIET") || ($scoring_method eq "PCT"))){
                if (exists $database->{MOTIFS}->[$j]->{POSITION_SCORING_MATRIX}) {
                    my $sum+=$_ foreach @{$database->{MOTIFS}->[$j]->{POSITION_SCORING_MATRIX}->[0]};
                    $sum = $sum + ($sum)**0.5;
                    push @counts,$sum;
                }
                elsif ((!exists $database->{MOTIFS}->[$j]->{POSITION_SCORING_MATRIX}) && ($arg{counts}==$defaults{counts})) {
                    die "Counts missing for motifs created in NestedMica.  $score_method requires counts.\n   Please pass parameters:  counts => x, where x is number of counts\n ";
                }
                else {
                    push @counts, $arg{counts};
                }
            }



            #If aligment method is NM then pass function the PPM's for both matrices and scoring method, along with counts(only needed if required by function)
            #Also, passes background frequencies.
            if ($alignment_method eq "global"){
                    $score=_global_align("PPM1"=>\@p_ppm1, "PPM2"=>\@p_ppm2,"scoring"=>$scoring_method,"background"=>\@background,"counts"=>\@counts);
            }

            elsif ($alignment_method eq "gapped"){
                    $score=_gapped_align("PPM1"=>\@p_ppm1, "PPM2"=>\@p_ppm2,"scoring"=>$scoring_method,"background"=>\@background,"counts"=>\@counts);
            }


            elsif ($alignment_method eq "trimmed"){
                    $score=_trimmed_align("PPM1"=>\@p_ppm1, "PPM2"=>\@p_ppm2,"scoring"=>$scoring_method,"background"=>\@background,"counts"=>\@counts, "trim"=>$trim);
            }

            else { die "Alignment Method not defined";}
            print "$query_label\t$database_label\t$score\n";
        }
    }

}


#############################################################################
#
#   score_motifset()
#
#   Description:  compare motifset and calculate distances
#
#   Work: call splitstree and weblogo from function
#
sub score_motifset{
    my %defaults=("score"=>"NM",
                  "align"=>"global",
                  "background"=>{'A'=>0.25, 'C'=>0.25, 'G'=>0.25, 'T'=>0.25},
                  "counts"=>undef,
                  "trim"=>4);
    my $self=shift;
    my %arg=(%defaults,@_);

    my $scoring_method=$arg{score};
    my $alignment_method=$arg{align};
    my $trim=$arg{trim};

    $arg{background}=($defaults{background} == $arg{background}) ? $arg{background}:  #if seq_freq then get single nucleotide frequency table
                    $arg{background}->{1};

    my @background;
    foreach my $letter (sort keys %{$arg{background}}){
        push @background, ${$arg{background}}{$letter};
    }

    my $align_method=   ($alignment_method eq "global") ? "global" :
                        ($alignment_method eq "gapped") ? "gapped"   :
                        ($alignment_method eq "trimmed") ? "trimmed"  :
                        die "Unsupported scoring method";

    my $score_method=   ($scoring_method eq "NM") ? "NestedMica" :
                        ($scoring_method eq "KLD") ? "Kullback-Leibler Distance" :
                        ($scoring_method eq "PCC") ? "Pearson's Correlation Coefficient" :
                        ($scoring_method eq "SW") ? "Sandelin/Wasserman" :
                        ($scoring_method eq "ALLR") ? "Average Log-likelihood Ratio" :
                        ($scoring_method eq "ALLR2") ? "Average Log-likelihood Ratio with limit of 2 imposed":
                        ($scoring_method eq "PCT") ?  "Pearson's Chi-Square Test" :
                        ($scoring_method eq "FIET") ? "Fisher-Irwin Exact Test" :
                        ($scoring_method eq "ED") ? "Euclidean Distance" :
                        die "Unsupported scoring method";

    my @motif_labels;
    foreach my $motif (@{$self->{MOTIFS}}){
        push @motif_labels,$motif->{ACCESSION};
    }

    my @score_table;
    for (my $i=0;$i<scalar(@{$self->{MOTIFS}});$i++){
        my @score_line;
        my $j=$i+1;
        my @p_ppm1=@{$self->{MOTIFS}->[$i]->{POSITION_PROBABILITY_MATRIX}};
        my @counts;

        #Find pseudocounted counts from scoring matrix, because these functions require counts.
        #First checks if the function requiring counts are called
        #Second checks to see if POSITION_SCORING_MATRIX is available for function.
        #If the PSM isn't available and counts weren't handed to function then die
        #If counts were handed to function then it uses these counts.
        if ((($scoring_method eq "ALLR") || ($scoring_method eq "ALLR2")||
             ($scoring_method eq "FIET") || ($scoring_method eq "PCT"))){
            if (exists $self->{MOTIFS}->[$i]->{POSITION_SCORING_MATRIX}) {
                my $sum+=$_ foreach @{$self->{MOTIFS}->[$i]->{POSITION_SCORING_MATRIX}->[0]};
                $sum = $sum + ($sum)**0.5;
                push @counts,$sum;
            }
            elsif ((!exists $self->{MOTIFS}->[$i]->{POSITION_SCORING_MATRIX}) && ($arg{counts}==$defaults{counts})) {
                die "Counts missing for motifs created in NestedMica.  $score_method requires counts.\n   Please pass parameters:  counts => x, where x is number of counts\n ";
            }
            else {
                push @counts, $arg{counts};
            }
        }

        for (my $j=0;$j<=$i;$j++){
            my @p_ppm2=@{$self->{MOTIFS}->[$j]->{POSITION_PROBABILITY_MATRIX}};

            #Find pseudocounted counts from scoring matrix, because these functions require counts.
            #First checks if the function requiring counts are called
            #Second checks to see if POSITION_SCORING_MATRIX is available for function.
            #If the PSM isn't available and counts weren't handed to function then die
            #If counts were handed to function then it uses these counts.
            if ((($scoring_method eq "ALLR") || ($scoring_method eq "ALLR2")||
                 ($scoring_method eq "FIET") || ($scoring_method eq "PCT"))){
                if (exists $self->{MOTIFS}->[$i]->{POSITION_SCORING_MATRIX}) {
                    my $sum+=$_ foreach @{$self->{MOTIFS}->[$i]->{POSITION_SCORING_MATRIX}->[0]};
                    $sum = $sum + ($sum)**0.5;
                    push @counts,$sum;
                }
                elsif ((!exists $self->{MOTIFS}->[$i]->{POSITION_SCORING_MATRIX}) && ($arg{counts}==$defaults{counts})) {
                    die "Counts missing for motifs created in NestedMica.  $score_method requires counts.\n   Please pass parameters:  counts => x, where x is number of counts\n ";
                }
                else {
                    push @counts, $arg{counts};
                }
            }



            #If aligment method is NM then pass function the PPM's for both matrices and scoring method, along with counts(only needed if required by function)
            #Also, passes background frequencies.
            if ($alignment_method eq "global"){
                    my $score=_global_align("PPM1"=>\@p_ppm1, "PPM2"=>\@p_ppm2,"scoring"=>$scoring_method,"background"=>\@background,"counts"=>\@counts);
                    push @score_line,$score;
            }

            elsif ($alignment_method eq "gapped"){
                    my $score=_gapped_align("PPM1"=>\@p_ppm1, "PPM2"=>\@p_ppm2,"scoring"=>$scoring_method,"background"=>\@background,"counts"=>\@counts);
                    push @score_line,$score;
            }


            elsif ($alignment_method eq "trimmed"){
                    my $score=_trimmed_align("PPM1"=>\@p_ppm1, "PPM2"=>\@p_ppm2,"scoring"=>$scoring_method,"background"=>\@background,"counts"=>\@counts, "trim"=>$trim);
                    push @score_line,$score;
            }

            else { die "Alignment Method not defined";}
        }
        push @score_table,\@score_line;
    }

    my %score=("ALIGN_METHOD" =>$alignment_method,
               "SCORING_METHOD"=>$scoring_method,
               "DISTANCE_MATRIX" => \@score_table,
               "LABEL" => \@motif_labels);
    $self->{DISTANCE_SCORE}=\%score;
    print "Distance_Score created using methods: \n Align: $alignment_method\n Score: $scoring_method\n";
    return $self;
}


###############################################################################
#
#   _trim_align
#
#   trimmed-local alignment compares
#
#
sub _trimmed_align{
   my %defaults=("PPM1"=>undef,
                 "PPM2"=>undef,
                 "scoring"=>"NM",
                 "background"=>undef,
                 "counts"=>undef,
                 "offset"=>undef,
                 "current_sum"=>undef,
                 "trim"=>undef);

    my %arg=(%defaults,@_);
    my $r_counts=$arg{counts};
    my $trim=$arg{trim};

    my @self1=@{$arg{PPM1}};
    my @self2=@{$arg{PPM2}};
    my $size1=scalar @self1;
    my $size2=scalar @self2;


    if (($size1<$trim)||($size2<$trim)){
        return "NA";
        #die "Motif is smaller than trimmed length";
    }

    #determine high information core compared to background
    my @KL;
    for (my $i=0;$i<$size1;$i++){
        my $KL=0;
        for (my $j=0;$j<4;$j++){
            $KL+=($self1[$i][$j]*log($self1[$i][$j]/$arg{background}->[$j]))/log(2);
        }
        push @KL, $KL;
    }

    my $position=0;
    my $top_sum=0;
    for (my $i=0;$i<=$size1-$trim;$i++){
        my $sum=0;
        for (my $j=$i;$j<$i+$trim;$j++){
            $sum+=$KL[$j];
        }
        if ($sum>$top_sum){
            $position=$i;
            $top_sum=$sum;
        }
    }

    #restrict self1 to core position
    @self1=@self1[$position..$position+$trim-1];

    #scan self2 with core position of $self1(foreach possible alignment that is trimmed to specified length)

    for (my $i=0;$i<$size2-$trim+1;$i++){
        my @ppm2=@self2[$i..$i+$trim-1];

        my @rev_self2;
        for (my $i=(scalar @ppm2)-1;$i>=0;$i--){
            my @reversed=reverse @{$ppm2[$i]};
            push @rev_self2, \@reversed;
        }

        my $new_sum;
        my $new_rev_sum;
        if ($arg{scoring} eq "NM") {
            $new_sum=_nm(\@self1,\@ppm2);   #calculate score
            $new_rev_sum=_nm(\@self1,\@ppm2);
        }
        elsif ($arg{scoring} eq "KLD") {
            $new_sum=_kld(\@self1,\@ppm2);   #calculate score
            $new_rev_sum=_kld(\@self1,\@rev_self2);
        }
        elsif ($arg{scoring} eq "ED"){
            $new_sum=_ed(\@self1,\@ppm2);   #calculate score
            $new_rev_sum=_ed(\@self1,\@rev_self2);
        }
        elsif ($arg{scoring} eq "PCC"){
            $new_sum=_pcc(\@self1,\@ppm2);   #calculate score
            $new_rev_sum=_pcc(\@self1,\@rev_self2);
        }
        elsif ($arg{scoring} eq "ALLR"){
            $new_sum=_allr(\@self1,\@ppm2, $r_counts,$arg{background});   #calculate score
            $new_rev_sum=_allr(\@self1,\@rev_self2, $r_counts,$arg{background});
        }
        elsif ($arg{scoring} eq "ALLR2"){
            $new_sum=_allr2(\@self1,\@ppm2, $r_counts,$arg{background});   #calculate score
            $new_rev_sum=_allr2(\@self1,\@rev_self2, $r_counts,$arg{background});
        }
        elsif ($arg{scoring} eq "PCT"){
            $new_sum=_pct(\@self1,\@ppm2, $r_counts);   #calculate score
            $new_rev_sum=_pct(\@self1,\@rev_self2, $r_counts);
        }
        elsif ($arg{scoring} eq "FIET"){
            $new_sum=_fiet(\@self1,\@ppm2, $r_counts);   #calculate score
            $new_rev_sum=_fiet(\@self1,\@rev_self2, $r_counts);
        }
        elsif ($arg{scoring} eq "SW"){
            $new_sum=_sw(\@self1,\@ppm2);   #calculate score
            $new_rev_sum=_sw(\@self1,\@rev_self2);
        }
        elsif ($arg{scoring} eq "ESS"){
            $new_sum=_ess(\@self1,\@ppm2);   #calculate score
            $new_rev_sum=_ess(\@self1,\@rev_self2);
        }
        else {
            die "Unknown Scoring Method";
        }

        ##Determine which score is closer to zero and return the closer one.
        $new_sum=(abs($new_rev_sum)<abs($new_sum)) ? $new_rev_sum : $new_sum;  #Return smaller sum of Reverse or Forward
        if (defined $arg{current_sum}){
            $arg{current_sum}= (abs($new_sum)<abs($arg{current_sum})) ? $new_sum : $arg{current_sum}; }
        else {$arg{current_sum}=$new_sum;}#if $new_sum is smaller save to $sum

    }

    return $arg{current_sum};
}

###############################################################################
#
#   _gapped_align
#
#   gapped-global alignment compares both with gaps
#
# TODO:  OPTIMIZE FUNCTION
sub _gapped_align{
   my %defaults=("PPM1"=>undef,
                 "PPM2"=>undef,
                 "scoring"=>"NM",
                 "background"=>undef,
                 "counts"=>undef,
                 "offset"=>undef,
                 "current_sum"=>undef);
    my %arg=(%defaults,@_);

    my @self1=@{$arg{PPM1}};
    my @self2=@{$arg{PPM2}};

    my $size1=scalar @self1;
    my $size2=scalar @self2;

    my $new_sum;

    for (my $i=0;$i<$size1-2;$i++){
        my @ppm1=@self1;
        my @ppm2=@self2;
        if ($i>0){
            splice(@ppm1,$i,1);
        }
        $new_sum=_global_align("PPM1"=>\@ppm1, "PPM2"=>\@ppm2,"scoring"=>$arg{scoring},"background"=>$arg{background},"counts"=>$arg{counts});

        if (defined $arg{current_sum}){
            $arg{current_sum}= (abs($new_sum)<abs($arg{current_sum})) ? $new_sum : $arg{current_sum}; }
        else {$arg{current_sum}=$new_sum;}#if $new_sum is smaller save to $sum
    }

    for (my $i=0;$i<$size2-2;$i++){
        my @ppm1=@self1;
        my @ppm2=@self2;
        if ($i>0){
            splice(@ppm2,$i,1);
        }
        $new_sum=_global_align("PPM1"=>\@ppm1, "PPM2"=>\@ppm2,"scoring"=>$arg{scoring},"background"=>$arg{background},"counts"=>$arg{counts});

        if (defined $arg{current_sum}){
            $arg{current_sum}= (abs($new_sum)<abs($arg{current_sum})) ? $new_sum : $arg{current_sum}; }
        else {$arg{current_sum}=$new_sum;}#if $new_sum is smaller save to $sum
    }

    return $arg{current_sum};
}

###############################################################################
#
#   _global_align
#   Nested Mica alignemnt scheme
#   (refer to http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=17238282&dopt=Abstract)
#
#

sub _global_align{
   my %defaults=("PPM1"=>undef,
                 "PPM2"=>undef,
                 "scoring"=>"NM",
                 "background"=>undef,
                 "counts"=>undef,
                 "offset"=>undef,
                 "current_sum"=>undef);
    my %arg=(%defaults,@_);
    my $r_counts=$arg{counts};

    my @self1=@{$arg{PPM1}};
    my @self2=@{$arg{PPM2}};
    if ((defined $arg{offset})&&($arg{offset}<(scalar(@{$arg{PPM2}})-1))){

        #Reverse compliment score
        my @rev_self2;
        for (my $i=(scalar @self2)-1;$i>=0;$i--){
            my @reversed=reverse @{$self2[$i]};
            push @rev_self2, \@reversed;
        }

        if ($arg{offset}<=0){
            for (my $i=$arg{offset};$i<0;$i++){
                unshift @self2, $arg{background};#Add background values to beginning of @self2 ->bigger than @self1
                unshift @rev_self2, $arg{background};}}
        elsif ($arg{offset}>0){
            for (my $i=$arg{offset};$i>0;$i--){
                unshift @self1,$arg{background};}}
        my $seq2offset=scalar(@self2)-scalar(@self1);  #Checks to see if @self1 is still longer than $self2
        if ($seq2offset>0){                         #if it is smaller then add some to @self1
            for (my $i=$seq2offset;$i>0;$i--){
                push @self1, $arg{background};}}
        elsif ($seq2offset<0){                      #else add it to @self2
            for (my $i=$seq2offset;$i<0;$i++){
                push @self2, $arg{background};
                push @rev_self2, $arg{background};}}

        my $new_sum;
        my $new_rev_sum;
        if ($arg{scoring} eq "NM") {
            $new_sum=_nm(\@self1,\@self2);   #calculate score
            $new_rev_sum=_nm(\@self1,\@rev_self2);
        }
        elsif ($arg{scoring} eq "KLD") {
            $new_sum=_kld(\@self1,\@self2);   #calculate score
            $new_rev_sum=_kld(\@self1,\@rev_self2);
        }
        elsif ($arg{scoring} eq "ED"){
            $new_sum=_ed(\@self1,\@self2);   #calculate score
            $new_rev_sum=_ed(\@self1,\@rev_self2);
        }
        elsif ($arg{scoring} eq "PCC"){
            $new_sum=_pcc(\@self1,\@self2);   #calculate score
            $new_rev_sum=_pcc(\@self1,\@rev_self2);
        }
        elsif ($arg{scoring} eq "ALLR"){
            $new_sum=_allr(\@self1,\@self2, $r_counts,$arg{background});   #calculate score
            $new_rev_sum=_allr(\@self1,\@rev_self2, $r_counts,$arg{background});
        }
        elsif ($arg{scoring} eq "ALLR2"){
            $new_sum=_allr2(\@self1,\@self2, $r_counts,$arg{background});   #calculate score
            $new_rev_sum=_allr2(\@self1,\@rev_self2, $r_counts,$arg{background});
        }
        elsif ($arg{scoring} eq "PCT"){
            $new_sum=_pct(\@self1,\@self2, $r_counts);   #calculate score
            $new_rev_sum=_pct(\@self1,\@rev_self2, $r_counts);
        }
        elsif ($arg{scoring} eq "FIET"){
            $new_sum=_fiet(\@self1,\@self2, $r_counts);   #calculate score
            $new_rev_sum=_fiet(\@self1,\@rev_self2, $r_counts);
        }
        elsif ($arg{scoring} eq "SW"){
            $new_sum=_sw(\@self1,\@self2);   #calculate score
            $new_rev_sum=_sw(\@self1,\@rev_self2);
        }
        elsif ($arg{scoring} eq "ESS"){
            $new_sum=_ess(\@self1,\@self2);   #calculate score
            $new_rev_sum=_ess(\@self1,\@rev_self2);
        }
        else {
            die "Unknown Scoring Method";
        }

        ##Determine which score is closer to zero and return the closer one.
        $new_sum=(abs($new_rev_sum)<abs($new_sum)) ? $new_rev_sum : $new_sum;  #Return smaller sum of Reverse or Forward
        if (defined $arg{current_sum}){
            $arg{current_sum}= (abs($new_sum)<abs($arg{current_sum})) ? $new_sum : $arg{current_sum}; }
        else {$arg{current_sum}=$new_sum;}#if $new_sum is smaller save to $sum

        $arg{offset}++;
        $arg{current_sum}=_global_align(%arg);}  #recursively call function for other offsets;}

    elsif (!defined $arg{offset}) {
        $arg{offset}=1-scalar(@self1);
        _global_align(%arg);}

    else {return $arg{current_sum};}
}

#Scoring Methods Below :
###############################################################################
#
#   _nm()
#
#   Description:  NestedMica scoring
#
#
#   Work:   Explain scoring and pros/cons
#
sub _nm{
    my ($r_array1,$r_array2)=@_;
    my $sum;
    for (my $i=0;$i<scalar(@{$r_array1});$i++){
        for (my $j=0;$j<4;$j++){
            $sum+= ((($r_array1->[$i]->[$j]-$r_array2->[$i]->[$j])**2)**(1.25));
        }
    }
    return $sum;
}
###############################################################################
#
#   _pcc()
#
#   Description:  Pearson Correlation Coefficient
#
#
#   Work:   Explain scoring and pros/cons
#
#
sub _pcc {
    my ($r_array1,$r_array2)=@_;
    my $sum;
    for (my $i=0;$i<scalar(@{$r_array1});$i++){
        my $PCC_denominator;
        my $PCC_numerator=0;
        for (my $j=0;$j<4;$j++){
            $PCC_numerator += ($r_array1->[$i]->[$j] - 0.25) * ($r_array2->[$i]->[$j] - 0.25);
            $PCC_denominator+=(($r_array1->[$i]->[$j] - 0.25)**2)*(($r_array2->[$i]->[$j] - 0.25)**2)
        }
        if ($PCC_denominator==0){
            next; }
        else {
            $sum+=($PCC_numerator/($PCC_denominator)**.5)-1; }
    }
    return $sum;
}

###############################################################################
#
#   _kld()
#
#   Description:  Kullback-Leibler Distance (Symmetrical distance)
#
#   Work:   Explain scoring and pros/cons
#
sub _kld{
    my ($r_array1,$r_array2)=@_;
    my $sum;
    for (my $i=0;$i<scalar(@{$r_array1});$i++){
        for (my $j=0;$j<4;$j++){
                $sum+=.5*($r_array1->[$i]->[$j]*((log($r_array1->[$i]->[$j]/$r_array2->[$i]->[$j])/log(2)))+
                      $r_array2->[$i]->[$j]*((log($r_array2->[$i]->[$j]/$r_array1->[$i]->[$j])/log(2))));
        }

    }
    return $sum;
}


###############################################################################
#
#   _allr()
#
#   Description:  Average Log-Likelihood
#
#   Work:   Explain scoring and pros/cons
#
#
sub _allr{
    my ($r_array1,$r_array2,$r_counts,$background)=@_;   #r_counts is a reference to matrix (total counts per position for each array);  (N,M)

    my $sum;
    for (my $i=0;$i<scalar(@{$r_array1});$i++){
        my $array1;
        my $array2;
        my $denominator;
        for (my $j=0;$j<4;$j++){
                $array1+=($r_array1->[$i]->[$j]*$r_counts->[0])*log($r_array1->[$i]->[$j]/$background->[$j]);
                $array2+=($r_array2->[$i]->[$j]*$r_counts->[1])*log($r_array2->[$i]->[$j]/$background->[$j]);
        }
        $denominator=($r_counts->[0]+$r_counts->[1]);
        $sum+=($array1+$array2)/$denominator;
    }
    return $sum;
}

###############################################################################
#
#   _allr2()
#
#   Description:  Average Log-Likelihood with a limit of -2 imposed
#
#   Work:   Explain scoring and pros/cons
#
#
sub _allr2 {
    my ($r_array1,$r_array2,$r_counts,$background)=@_;   #r_counts is a reference to matrix (total counts per position for each array);  (N,M)
    my $limited_sum;
    my $sum;
    for (my $i=0;$i<scalar(@{$r_array1});$i++){
        my $array1;
        my $array2;
        my $denominator;
        for (my $j=0;$j<4;$j++){
                $array1+=($r_array1->[$i]->[$j]*$r_counts->[0])*log($r_array1->[$i]->[$j]/$background->[$j]);
                $array2+=($r_array2->[$i]->[$j]*$r_counts->[1])*log($r_array2->[$i]->[$j]/$background->[$j]);
        }
        $denominator=($r_counts->[0]+$r_counts->[1]);
        $limited_sum=($array1+$array2)/$denominator;
        if (($sum+$limited_sum)<-2){
            $sum=-2;
        }
        else{
            $sum+=$limited_sum;
        }
    }
    return $sum;
}

###############################################################################
#
#   _sw()
#
#   Description:    Sandelin/Wasserman  (also refered to Sum of Squared distances)
#
#
#   Work:   Explain scoring and pros/cons
#
#
sub _sw{
    my ($r_array1,$r_array2)=@_;
    my $sum;
    for (my $i=0;$i<scalar(@{$r_array1});$i++){
        my $SW=0;
        for (my $j=0;$j<4;$j++){
                $SW+=($r_array1->[$i]->[$j] - $r_array2->[$i]->[$j])**2
        }
        #$sum+=2-$SW;
        $sum+=$SW;
    }
    return $sum;
}

###############################################################################
#
#   _pct()
#
#   Description:  Pearson Chi-squared test   (Array's handed to this function should be PSM)
#
#   Work:   Explain scoring and pros/cons
#
#
sub _pct {
    my ($r_array1,$r_array2,$r_counts)=@_;  #r_counts is a reference to matrix (total counts per position for each array);  (N,M)
    my $GAMMA=0.8862269255;
    my $pvalue=1;
    my $N=$r_counts->[1]+$r_counts->[0];  #Total counts for both arrays
    for (my $i=0;$i<scalar(@{$r_array1});$i++){
        my $score=0;
        #Contigency Table for column
        my @Nn=(($r_counts->[0]*$r_array1->[$i]->[0]+$r_array2->[$i]->[0]*$r_counts->[1]),
                ($r_counts->[0]*$r_array1->[$i]->[1]+$r_array2->[$i]->[1]*$r_counts->[1]),
                ($r_counts->[0]*$r_array1->[$i]->[2]+$r_array2->[$i]->[2]*$r_counts->[1]),
                ($r_counts->[0]*$r_array1->[$i]->[3]+$r_array2->[$i]->[3]*$r_counts->[1]));
        for (my $j=0;$j<4;$j++){
                my $N1j_exp=($r_counts->[0]*$Nn[$j])/$N;
                my $N2j_exp=($r_counts->[1]*$Nn[$j])/$N;
                if (($N1j_exp<5) | ($N2j_exp<5)){
                    die " Pearson's Chi-squared requires Expected values to be greater or equal to 5.\n";
                }

                $score+=((($r_array1->[$i]->[$j]*$r_counts->[0])-$N1j_exp)**2)/$N1j_exp;
                $score+=((($r_array2->[$i]->[$j]*$r_counts->[1])-$N2j_exp)**2)/$N1j_exp;
        }



        #Calculate  Incomplete Gamma from approximation (see http://www.rskey.org/gamma.htm)
        my $fact_value=1.32934038818;  #Factorial Value initialized at a=1.5
        my $exp=2.71828182846;  #eulers constant
        my $x=$score/2;  # score/2 for the input int incomplete gamma for Chi CDf
        my $sum=0;
        for (my $n=0;$n<=50;$n++){
            $sum+=($x**$n)/($fact_value);
            $fact_value=$fact_value*($fact_value+1);
        }
        my $IGAMMA=($x**1.5)*($exp**-$x)*$sum;


        #Calculate P-value for Chi score
        my $pvalue*=1-($IGAMMA/$GAMMA);
    }
    return $pvalue;
}

###############################################################################
#
#   _fiet()
#
#   Description:    Fisher-Irwin Exact Test  (Array handed to this function should be PSM function)
#                   Score is the log base 2 of probability returned from Fisher-Irwin exact test.
#
#   Work:   Explain scoring and pros/cons
#
#
# Warning Fisher-Irwin Exact test not optimized algorithm  => takes a long time
#
sub _fiet{
    my ($r_array1,$r_array2,$r_counts)=@_;  #r_counts is a reference to matrix (total counts per position for each array);  (N,M)
    my $N=$r_counts->[1]+$r_counts->[0];  #Total counts for both arrays
    my $score=0;
    for (my $i=0;$i<scalar(@{$r_array1});$i++){
        my @row1=($r_array1->[$i]->[0]*$r_counts->[0],
                  $r_array1->[$i]->[1]*$r_counts->[0],
                  $r_array1->[$i]->[2]*$r_counts->[0],
                  $r_array1->[$i]->[3]*$r_counts->[0]);
        my @row2=($r_array2->[$i]->[0]*$r_counts->[1],
                  $r_array2->[$i]->[1]*$r_counts->[1],
                  $r_array2->[$i]->[2]*$r_counts->[1],
                  $r_array2->[$i]->[3]*$r_counts->[1]);
        $score+=log(_fisher(\@row1,\@row2))/log(2);
    }
    return $score;
}


#Fisher-Irwin Exact Test
#Input is two columns of scores
#Fisher-Irwin Exact test for 2x4 contigency table
#also known as Freeman-Halton  extension of the Fisher exact probability test for a two-rows by four-columns
#Directly translated from http://faculty.vassar.edu/lowry/fisher2x4.html
#If expected number is larger than 5 it is better to use the Chi-Squared test

#Fisher exact test is used where total size of the data set are smaller than 120.
#Results are Non-directional two-tailed probabilities

#Fisher exact test takes counts, therefore they must be rounded up or down.
sub _fisher{
    my ($r_row1, $r_row2)=@_;


    #Round function:  if above x.5 then rounds up if less than x.5 then rounds down.
    my $round=sub {  #
        my($number) = shift;
        return int($number + .5);
        };

    my @tabx=(@{$r_row1},@{$r_row2});

    @tabx=map &$round($_),@tabx;  #Round values in table.

    my ($n,$r1,$r2,$c1,$c2,$c3,$c4,@taby,@tabj,@expc,@cscomp,$holdx1,$holdx2);
    my $totprob;
    my $totprob2;
    @tabj=@tabx;

    $r1=$tabx[0]+$tabx[1]+$tabx[2]+$tabx[3];
    $r2=$tabx[4]+$tabx[5]+$tabx[6]+$tabx[7];

    $c1=$tabx[0]+$tabx[4];
    $c2=$tabx[1]+$tabx[5];
    $c3=$tabx[2]+$tabx[6];
    $c4=$tabx[3]+$tabx[7];

    $n=$r1+$r2;
    my $eq=0;
    if (($c1==$c2) && ($c2==$c3)){$eq+=1;}
    if ($c3==$c4){$eq+=1;}

    if ($eq>1){
        for (my $i=0;$i<8;$i++){
            if (abs($tabj[1]-$tabj[5])>abs($tabj[0]-$tabj[4])){
                $holdx1=$tabj[0];
                $holdx2=$tabj[4];
                $tabj[0]=$tabj[1];
                $tabj[4]=$tabj[5];
                $tabj[1]=$holdx1;
                $tabj[5]=$holdx2;
            }
            elsif (abs($tabj[2]-$tabj[6])>abs($tabj[1]-$tabj[5])){
                $holdx1=$tabj[1];
                $holdx2=$tabj[5];
                $tabj[1]=$tabj[2];
                $tabj[5]=$tabj[6];
                $tabj[2]=$holdx1;
                $tabj[6]=$holdx2;
            }
            elsif (abs($tabj[3]-$tabj[7])>abs($tabj[2]-$tabj[6])){
                $holdx1=$tabj[2];
                $holdx2=$tabj[6];
                $tabj[2]=$tabj[3];
                $tabj[6]=$tabj[7];
                $tabj[3]=$holdx1;
                $tabj[7]=$holdx2;
            }
        }
    }
    else {
        for(my $i=0;$i<8;$i++){
            if (($tabj[1]+$tabj[5])>($tabj[0]+$tabj[4])){
                $holdx1=$tabj[0];
                $holdx2=$tabj[4];
                $tabj[0]=$tabj[1];
                $tabj[4]=$tabj[5];
                $tabj[1]=$holdx1;
                $tabj[5]=$holdx2;
            }
            elsif (($tabj[2]+$tabj[6])>($tabj[1]+$tabj[5])){
                $holdx1=$tabj[1];
                $holdx2=$tabj[5];
                $tabj[1]=$tabj[2];
                $tabj[5]=$tabj[6];
                $tabj[2]=$holdx1;
                $tabj[6]=$holdx2;
            }
            elsif (($tabj[3]+$tabj[7])>($tabj[2]+$tabj[6])){
                $holdx1=$tabj[2];
                $holdx2=$tabj[6];
                $tabj[2]=$tabj[3];
                $tabj[6]=$tabj[7];
                $tabj[3]=$holdx1;
                $tabj[7]=$holdx2;
            }
        }
    }

    $r1 = $tabx[0]+$tabx[1]+$tabx[2]+$tabx[3];
    $r2 = $tabx[4]+$tabx[5]+$tabx[6]+$tabx[7];

    $c1 = $tabx[0]+$tabx[4];
    $c2 = $tabx[1]+$tabx[5];
    $c3 = $tabx[2]+$tabx[6];
    $c4 = $tabx[3]+$tabx[7];
    $n = $r1+$r2;


    my $count = 0;
    my $tot = 0;

    if ($n<=120){
        @tabx=@tabj;
        $r1 = $tabx[0]+$tabx[1]+$tabx[2]+$tabx[3];
        $r2 = $tabx[4]+$tabx[5]+$tabx[6]+$tabx[7];
        $c1 = $tabx[0]+$tabx[4];
        $c2 = $tabx[1]+$tabx[5];
        $c3 = $tabx[2]+$tabx[6];
        $c4 = $tabx[3]+$tabx[7];
        $n = $r1+$r2;

        my ($factn,$facta,$factb,$factc,$factd,$facte,$factf,$factg,$facth,$factr1,$factr2,$factc1,$factc2,$factc3,$factc4)=
            (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

        @taby=@tabx;


        for(my $i=1;$i<=$n;$i++){$factn+=log($i)};

        for(my $i=1;$i<=$tabx[0];$i++){$facta+=log($i)};
        for(my $i=1;$i<=$tabx[1];$i++){$factb+=log($i)};
        for(my $i=1;$i<=$tabx[2];$i++){$factc+=log($i)};
        for(my $i=1;$i<=$tabx[3];$i++){$factd+=log($i)};
        for(my $i=1;$i<=$tabx[4];$i++){$facte+=log($i)};
        for(my $i=1;$i<=$tabx[5];$i++){$factf+=log($i)};
        for(my $i=1;$i<=$tabx[6];$i++){$factg+=log($i)};
        for(my $i=1;$i<=$tabx[7];$i++){$facth+=log($i)};
        for(my $i=1;$i<=$r1;$i++){$factr1+=log($i)};
        for(my $i=1;$i<=$r2;$i++){$factr2+=log($i)};
        for(my $i=1;$i<=$c1;$i++){$factc1+=log($i)};
        for(my $i=1;$i<=$c2;$i++){$factc2+=log($i)};
        for(my $i=1;$i<=$c3;$i++){$factc3+=log($i)};
        for(my $i=1;$i<=$c4;$i++){$factc4+=log($i)};

        my $num = $factr1+$factr2+$factc1+$factc2+$factc3+$factc4;
        my $denom = $factn+$facta+$factb+$factc+$factd+$facte+$factf+$factg+$facth;
        my $pex = exp($num-$denom);
        $totprob = 0;
        $totprob2 = 0;
        my $temp = 0;

        for(my $k=0; $k<($c1+1);$k++){
            $taby[0] = $k;
            $taby[4] = $c1-$k;
            for(my $j=0; $j<($c2+1);$j++){
                $taby[1] = $j;
                $taby[5] = $c2-$j;
                for(my $m=0; $m<($c3+1);$m++){
                    $taby[2] = $m;
                    $taby[6] = $c3-$m;
                    $taby[3] = $r1-$taby[0]-$taby[1]-$taby[2];
                    $taby[7] = $r2-$taby[4]-$taby[5]-$taby[6];
                    $tot = 0;
                    for(my $l=0; $l<8;$l++){$tot+=abs($taby[$l])}
                    if($tot==$n){
                        $count++;
                        $facta=0;
                        $factb=0;
                        $factc=0;
                        $factd=0;
                        $facte=0;
                        $factf=0;
                        $factg=0;
                        $facth=0;
                        for(my $i=1; $i<=$taby[0]; $i++){$facta+=log($i)};
                        for(my $i=1; $i<=$taby[1]; $i++){$factb+=log($i)};
                        for(my $i=1; $i<=$taby[2]; $i++){$factc+=log($i)};
                        for(my $i=1; $i<=$taby[3]; $i++){$factd+=log($i)};
                        for(my $i=1; $i<=$taby[4]; $i++){$facte+=log($i)};
                        for(my $i=1; $i<=$taby[5]; $i++){$factf+=log($i)};
                        for(my $i=1; $i<=$taby[6]; $i++){$factg+=log($i)};
                        for(my $i=1; $i<=$taby[7]; $i++){$facth+=log($i)};
                        $denom = $factn+$facta+$factb+$factc+$factd+$facte+$factf+$factg+$facth;
                        $temp = exp($num-$denom);
                        if ($temp<=$pex){$totprob+=$temp}
                        if ($temp<$pex){$totprob2+=$temp}
                    }
                }
            }
        }

        $totprob2 = $totprob2+$pex;

        if ($totprob>=1){$totprob="1.0"}
        if ($totprob2>=1){$totprob2="1.0"}


        #print "Total PROB 1: $totprob\n";  #totprob= the probability of the observed array of cell frequencies plus the sum of the probabilities of all other cell-frequency arrays (such as would be consistent with the observed marginal totals) that are equal to or smaller than the probability of the observed array.
        #print "Total PROB 2: $totprob2\n";  #totprob2 = the probability of the observed array of cell frequencies plus the sum of the probabilities of all other cell-frequency arrays (such as would be consistent with the observed marginal totals) that are smaller than the probability of the observed array.
        #print $count;
    }
    if ((1-$totprob)<=0.000001){
        $totprob=1;
    }
 return $totprob;
}

###############################################################################
#
#   _ess()
#
#   Description:    Evolutionary Substitution Score
#
#   Work:   Explain scoring and pros/cons
#
#
sub _ess{
    my ($r_array1,$r_array2)=@_;
    my $sum;
    for (my $i=0;$i<scalar(@{$r_array1});$i++){
        for (my $j=0;$j<4;$j++){
                #function goes here
        }

    }
    return $sum;
}


###############################################################################
#
#   _ed()
#
#   Description:    Euclidian Distance
#                   Sum Euclidian distance for each column
#
#   Work:   Explain scoring and pros/cons
#
#
sub _ed{
    my ($r_array1,$r_array2)=@_;
    my $sum;
    for (my $i=0;$i<scalar(@{$r_array1});$i++){
        my $ed;
        for (my $j=0;$j<4;$j++){
                $ed=($r_array1->[$i]->[$j] - $r_array2->[$i]->[$j])**2;
        }
        $sum+=($ed)**.5;
    }
    return $sum;
}


################################################################################
#
#   splitstree
#
#   Description:  creates nexus file.
#
#   Work:  Setup file to be handled from score_motifset function rather than
#          accessed separately.  Depreciate directory creation
#
#
sub splitstree{
    my ($self,$filename)=@_;

    my $size=scalar(@{$self->{MOTIFS}});  #Only used to deterimine if set is going to be to big for SplitsTree
    if ($size>=550){
        print "You may experience problems with SplitsTree.  It doesn't handle this large of a dataset.\n";
    }

    my $Distance_Matrix=$self->{DISTANCE_SCORE}->{DISTANCE_MATRIX};
    my $Labels=$self->{DISTANCE_SCORE}->{LABEL};
    my $score_method=$self->{DISTANCE_SCORE}->{SCORING_METHOD};
    my $num=scalar @$Labels;

    open NEXUS ,">$filename.nex"  || die "Cannot open file for writing $!";
    print "FILENAME: $filename.nex\n";
    print NEXUS  "#NEXUS\nBEGIN taxa;\n";
    print NEXUS "\tDIMENSIONS ntax=" . $num . ";\n";
    print NEXUS "TAXLABELS\n";
    for (my $i=1;$i<=$num;$i++){
        my $n=$i-1;
        print NEXUS  "[$i]\t$Labels->[$n]\n";
    }
    print NEXUS  ";\nEND;\nBEGIN distances;\n";
    print NEXUS  "\tDIMENSIONS ntax=" . scalar @$Labels . ";\n";
    print NEXUS  "\tFORMAT\n\t\ttriangle=LOWER\n\t\tdiagonal\n\t\tlabels\n\t;\n";
    print NEXUS  "\tMATRIX\n";
    for (my $i=0;$i<$num;$i++){
        print NEXUS "\t$Labels->[$i]\t";
        my $line=join("\t",@{$Distance_Matrix->[$i]});
        print NEXUS  "$line\n";
    }
    print NEXUS  "\t;\n";
    print NEXUS  "END;\n";
    $|=1;
    close NEXUS;
    END:
}


################################
#
#   draw_motif
#
#   Description:  creates seq_logo summary HTML page.
#
sub _draw_motif{
    my @color=("green","orange","blue","red");
    my @letters=("A","C","G","T");
    my %defaults=("ppm"=>undef,
                  "filename"=>undef,
                  "directory"=>undef,
                  "template"=>undef,
                  "type"=>"pictogram",
                  "color"=>\@color,
                  "letter"=>\@letters);

    my %arg=(%defaults,@_);
    my $ppm=$arg{ppm};
    my $size=scalar @$ppm;  #size of motif
    my $directory=(defined $arg{directory}) ? $arg{directory}: getcwd();
    my $filename=(defined $arg{filename}) ? $arg{filename}:  die "filename not defined";
    my $rev_filename="Rev" . $filename;  #set filename for reverse complement motif

    # If template is undefined then import it.   If it is defined then skip import
    # This cuts down on imports if recursively calling function to draw.
    if (!defined $arg{template}) {
        $arg{template}=main::svg_template();
    }


    my $graph_length=70+(80*$size);

    my @forward_svg=@{$arg{template}}[0..2];
    $forward_svg[2]=~s/\$length/$graph_length/;
    my @reverse_svg;


    my @x_axis_svg=@{$arg{template}}[4..16];
    my $adjusted_graph_length=$graph_length-30;
    $x_axis_svg[12]=~s/\$graph_length/$adjusted_graph_length/;
    $x_axis_svg[1]=~s/\$graphlength/$adjusted_graph_length/;

    if ($arg{type} eq "pictogram"){
        @reverse_svg=@forward_svg;
        my $horizontal_position=70;

        # set heights for pictogram
        for (my $i=0;$i<$size;$i++){

            my $character=$i+1;
            my $for_vertical_position=118;
            my $rev_vertical_position=118;


            my $axis=${$arg{template}}[18];
            my $position=$i+1;
            my $position_segment=70+($i*80);
            $axis=~s/\$position_segment/$position_segment/;
            $axis=~s/\$position/$position/;
            push @x_axis_svg, $axis;

            for (my $j=0;$j<=3;$j++){
                my $color=${$arg{color}}[$j];
                my $letter=$arg{letter}->[$j];
                if ($ppm->[$i]->[$j]>=0.01){
                    my $for_line=${$arg{template}}[3];
                    $for_line=~s/\$transformation/$ppm->[$i]->[$j]/;
                    $for_line=~s/\$horizontal_position/$horizontal_position/;
                    $for_line=~s/\$vertical_position/$for_vertical_position/;
                    $for_line=~s/\$letter_color/$color/;
                    $for_line=~s/\$letter/$letter/;
                    push @forward_svg, $for_line;
                    $for_vertical_position-=($ppm->[$i]->[$j]*88);
                }

                if ($ppm->[$size-1-$i]->[3-$j]>=0.01){
                    my $rev_line=${$arg{template}}[3];
                    $rev_line=~s/\$transformation/$ppm->[$size-1-$i]->[3-$j]/;
                    $rev_line=~s/\$horizontal_position/$horizontal_position/;
                    $rev_line=~s/\$vertical_position/$rev_vertical_position/;
                    $rev_line=~s/\$letter_color/$color/;
                    $rev_line=~s/\$letter/$letter/;
                    push @reverse_svg, $rev_line;
                    $rev_vertical_position-=($ppm->[$size-1-$i]->[3-$j]*88);
                }
            }
            $ horizontal_position+=80;
        }
    }

    elsif ($arg{type} eq "logo") {
        #change template entries for logo
        $x_axis_svg[9]=~s/Percentage/bits/;
        $x_axis_svg[9]=~s/110/75/;
        $x_axis_svg[8]=~s/100/2/;
        $x_axis_svg[8]=~s/6/18/;
        $x_axis_svg[7]=~s/50/1/;
        $x_axis_svg[7]=~s/12/18/;
        @reverse_svg=@forward_svg;
        my $horizontal_position=70;

        #calculate height of each base using Shannon's Information Content
        # H=-
        # set heights for logo
        for (my $i=0;$i<$size;$i++){

            my $character=$i+1;
            my $for_vertical_position=118;
            my $rev_vertical_position=118;


            my $axis=${$arg{template}}[18];
            my $position=$i+1;
            my $position_segment=70+($i*80);
            $axis=~s/\$position_segment/$position_segment/;
            $axis=~s/\$position/$position/;
            push @x_axis_svg, $axis;


            #Forward motif entropy
            my $for_h=0;  #Shannon's Information Content
            my $rev_h=0;
            for (my $j=0;$j<=3;$j++){
                $for_h+=-1*$ppm->[$i]->[$j]*(log($ppm->[$i]->[$j])/log(2));
                $rev_h+=-1*$ppm->[$size-1-$i]->[3-$j]*(log($ppm->[$size-1-$i]->[3-$j])/log(2));
            }

            #Total height of all letters
            my $for_total_height=2-$for_h;
            my $rev_total_height=2-$rev_h;

			# will store heights in order to put them in size order so biggest nucleotide is always on top
			my %forward_heights;
			my %reverse_heights;
			
			for (my $j=0;$j<=3;$j++){
	                if ($ppm->[$i]->[$j]>=0.01){
	                    my $for_line=${$arg{template}}[3];
	                    my $height=$ppm->[$i]->[$j]*$for_total_height/2;
						$forward_heights{$j}=$height;
	                }

	                if ($ppm->[$size-1-$i]->[3-$j]>=0.01){
	                    my $height=$ppm->[$size-1-$i]->[3-$j]*$rev_total_height/2;
						$reverse_heights{$j}=$height;
	                }
	         }			
				
	
	
			# loop through each nucleotide, starting with the lowest height 
			# first do forward motif...
			
			foreach my $j (sort {$forward_heights{$a} <=> $forward_heights{$b}}(keys %forward_heights)){
                my $color=${$arg{color}}[$j];
                my $letter=$arg{letter}->[$j];
                if ($ppm->[$i]->[$j]>=0.01){
                    my $for_line=${$arg{template}}[3];
                    my $height=$ppm->[$i]->[$j]*$for_total_height/2;
                    $for_line=~s/\$transformation/$height/;
                    $for_line=~s/\$horizontal_position/$horizontal_position/;
                    $for_line=~s/\$vertical_position/$for_vertical_position/;
                    $for_line=~s/\$letter_color/$color/;
                    $for_line=~s/\$letter/$letter/;
                    push @forward_svg, $for_line;
                    $for_vertical_position-=($height*88);
                }
			}
			
			#...then reverse strand motif
			foreach my $j (sort {$reverse_heights{$a} <=> $reverse_heights{$b}}(keys %reverse_heights)){
				my $color=${$arg{color}}[$j];
			    my $letter=$arg{letter}->[$j];
                if ($ppm->[$size-1-$i]->[3-$j]>=0.01){
                    my $rev_line=${$arg{template}}[3];
                    my $height=$ppm->[$size-1-$i]->[3-$j]*$rev_total_height/2;
                    $rev_line=~s/\$transformation/$height/;
                    $rev_line=~s/\$horizontal_position/$horizontal_position/;
                    $rev_line=~s/\$vertical_position/$rev_vertical_position/;
                    $rev_line=~s/\$letter_color/$color/;
                    $rev_line=~s/\$letter/$letter/;
                    push @reverse_svg, $rev_line;
                    $rev_vertical_position-=($height*88);
                }
            }
            $ horizontal_position+=80;
        }
    }
    else {
        die "$arg{type}: Incorrect type\n Possible types: pictogram, logo";
    }



    open (FSVG, "> $directory/$filename.svg");
    open (RSVG, "> $directory/$rev_filename.svg");
    print FSVG @forward_svg;
    print FSVG @x_axis_svg;
    print FSVG ${$arg{template}}[19];
    print RSVG @reverse_svg;
    print RSVG @x_axis_svg;
    print RSVG ${$arg{template}}[19];

    $|=1;
    close FSVG;
    close RSVG;

}


###############################
#
#   pictogram
#
#   Description:  creates seq_logo summary HTML page.
#
#   Work:  Modify files to handle different directory structure.
#           Work SVG file from splits to contain seq_logo files
#           Include option to do information content instead of probability
sub pictogram{
    my ($self,$job_directory)=@_;
    if (!defined $job_directory){
        die "No Directory for pictograms" }
    else {
    	mkdir("$job_directory",0755);
    }

    for(my $i=0; $i<scalar(@{$self->{MOTIFS}});$i++){

        my %defaults=(  "ppm"=>$self->{MOTIFS}->[$i]->{POSITION_PROBABILITY_MATRIX},
                        "filename"=>"$self->{MOTIFS}->[$i]->{ACCESSION}",
                        "directory"=>$job_directory,
                        "template"=>main::svg_template(),
                        "type"=>"pictogram");
        _draw_motif(%defaults);
    }

}

###############################
#
#   logo
#
#   Description:  creates seq_logo summary HTML page.
#
#   Work:  Modify files to handle different directory structure.
#           Work SVG file from splits to contain seq_logo files
#           Include option to do information content instead of probability
sub logo{
    my ($self,$job_directory)=@_;
    if (!defined $job_directory){
        die "No Directory for pictograms" }
    else {
    	mkdir("$job_directory",0755);
    }

    for(my $i=0; $i<scalar(@{$self->{MOTIFS}});$i++){

        my %defaults=(  "ppm"=>$self->{MOTIFS}->[$i]->{POSITION_PROBABILITY_MATRIX},
                        "filename"=>"$self->{MOTIFS}->[$i]->{ACCESSION}",
                        "directory"=>$job_directory,
                        "template"=>main::svg_template(),
                        "type"=>"logo");
        _draw_motif(%defaults);
    }

}


################################
#
#   html_summary
#
#   Description:  creates seq_logo summary HTML page.
#
#   Work:  Modify files to handle different directory structure.
#           Work SVG file from splits to contain seq_logo files
#           Include option to do information content instead of probability
sub html_summary{
    my ($self,$directory,$type)=@_;
    if (!defined $directory){
        die "No Directory for pictograms" }
    else {
    	mkdir("$directory",0755);
    	mkdir("$directory/LOGO",0755);
    }

    $type= (defined $type) ? $type : "pictogram";


    if ($type eq "pictogram"){
        pictogram($self,"$directory/LOGO");
    }
    else {
        logo($self,"$directory/LOGO");
    }

    my $HTML_iterator=0;
    #Creation of HTML Page for Job
    my @index_HTML=@{main::html_index()};
    my @summary_HTML_head=@{main::html_summary()};
    my @summary_HTML_body;
    my @menu_HTML=@{main::html_menu()};
    my @menu_HTML_links;



    foreach my $motif (@{$self->{MOTIFS}}){
        my $filename=$motif->{ACCESSION};
        my $rev_filename="Rev" . $filename;

        push @summary_HTML_body, ("\<h2\>\<a name\=\"$filename\"\>$filename\<\/a\>\<\/h2\>\n",
                          "\<table align\=\"center\"\>\n",
                          "\<tr\>\n",
                          "\<th\>Forward\<\/th\>\n",
                          "\<th\>Reverse\<\/th\>\n",
                          "\<\/tr\>\n",
                          "\<tr\>\n",
                          "\<td\>\<img src\=\"LOGO\/$filename.svg\" width\=\"306\" height\=\"84\" align\=\"middle\"\>\<\/td\>\n",
                          "\<td\>\<img src\=\"LOGO\/$rev_filename.svg\" width\=\"306\" height\=\"84\" align\=\"middle\"\>\<\/td\>\n",
                          "\<\/tr\>\n",
                          "\<\/table\>\n",
                          "\<br\>\n",
                          "\<hr\>");
        push @menu_HTML_links, "\<OPTION VALUE\=\"Summary.html\#$filename\"\>$filename\n";
        unlink "$directory/LOGO/$rev_filename.eps";
        unlink "$directory/LOGO/$filename.eps";
    }

    #Output Summary.html
    push @summary_HTML_body, ("\<\/body\>\n","\<\/html\>\n");
    open SUM_HTML, "> $directory/Summary.html";
    print SUM_HTML join("\n",@summary_HTML_head);
    print SUM_HTML join("\n",@summary_HTML_body);
    $|=1;
    close (SUM_HTML);
    #Output Index.html
    open INDEX, "> $directory/index.html";
    print INDEX join("\n",@index_HTML);
     $|=1;
    close (INDEX);

    #Output Menu's.html
    for (my $i=1; $i<6;$i++){
        my $name="Menu". $i . "\.html";
        my $Target="Motif$i";
        my @Target=@menu_HTML[0..29];
        $Target[27]=~s/\{\$\S+\}/$Target\n/;
        open MENU, "> $directory/$name";
        print MENU join("\n",@Target);
        print MENU join("\n",@menu_HTML_links);
        print MENU join("\n",@menu_HTML[30..34]);
        $|=1;
        close (MENU);
    }
    $|=1
}



############################################################################
#  Transfac_import
#
#  Description:    Imports transfac matrix file into Motifset object.   Parses each
#                  line of transfac file according to line identifier, which is the
#                  first two letters or numbers of a line.  PWM is created from the
#                  Postion scoring matrix (PSM) by pseudocounting the PSM and then
#                  calculating the log-likelihood ration according to a background
#                  frequency of 0.25 for each nucleotide.  Position
#                  Probability Matrix (PPM) is created from the PSM without pseudo-
#                  counting.
#
#  usage:  $Motifset->Transfac_import("Transfac File")
#
#  Example $Motifset->Transfac_import("Matrix.dat")
#      where Matrix.dat is the Transfac Matrix database.
#
#
#  Minimal matrix file information:
#  AC  M00001
#  XX
#  P0      A      C      G      T
#  01      1      2      2      0
#  02      2      1      2      0
#  03      3      0      1      1
#  04      0      5      0      0
#  05      5      0      0      0
#  06      0      0      4      1
#  07      0      1      4      0
#  08      0      0      0      5
#  09      0      0      5      0
#  10      0      1      2      2
#  11      0      2      0      3
#  12      1      0      3      1
#  XX
#  //
#
# End of matrix record is denoted by double backslash //
# Motif accession becomes TFM00001 ; denoting the matrix is in the Transfac format
#
#  You can append matrix files together
#
# Object Line identifiers:
#      AC :  Accession
#      ID :  Identifier
#      DT :  Date and Author
#      NA :  Name of binding site
#      DE :  Short Factor Description
#      BF :  List of linked factor entrees
#      PO :  Position Scoring Matrix
#      BA :  Statistical Basis
#      BS :  Binding Sites underlying Matrix
#      CC :  Comments
#      RX :  Medline ID
#      RN :  Reference #
#      RA :  Reference Authors
#      RT :  Reference Title
#      RL :  Reference Data
#      XX : delimiter between lines  (see above example)
#      // : end of Motif (see above example)
#
#
#
#
# PPM calculated from PSM
# PWM calculated from PSM (pseudocounted) and then calculating log-likelihood ratio
#
sub transfac_import {
    my ($self,$filename)=@_;

    #Delete any existing Distance Score from Motifset
    #b/c importing will change the motifset so score won't reflect current motifset
    if (exists $self->{Distance_Score}){
        delete $self->{Distance_Score};
    }

    my %transfac_header=(   "ID" => "IDENTIFIER",
                            "DT" => "DATE_AUTHOR_UPDATE_DATES",
                            "CO" => "COPYRIGHT",
                            "NA" => "NAME_OF_BINDING_SITE",
                            "DE" => "SHORT_FACTOR_DESCRIPTION",
                            "BF" => "LIST_OF_LINKED_FACTOR_ENTRIES",
                            "BA" => "STATISTICAL_BASIS",
                            "BS" => "BINDING_SITES_UNDERLYING_MATRIX",
                            "CC" => "COMMENTS",
                            "RX" => "MEDLINE_ID",
                            "RN" => "REFERENCE_#",
                            "RA" => "REFERENCE_AUTHORS",
                            "RT" => "REFERENCE_TITLE",
                            "RL" => "REFERENCE_DATA",
                            "MS" => "MOTIF_SEQUENCES",
                            "SE" => "SEQUENCES",
                            "MM" => "MOTIF_METHOD",
                            "ML" => "MOTIF_LOCATIONS",
                            "FI" => "FILE",
                            "CL" => "COMMAND_LINE",
                            "DF" => "DATAFILE",
                            "AL" => "ALPHABET",
                            "RE" => "REGULAR_EXPRESSION",
                            "CS" => "CONSENSUS_SEQUENCE",
                            "VV" => "VERSION",
                            "TH" => "THRESHOLD",
                            "JO" => "MOTIF_JOB",
                            "OR" => "ORGANISM",
                            "AS" => "OLD_ACCESSION");

    open (TRANSFAC, "< $filename") || die "Couldn't Open $filename:$!\n;";


    my $version;    #contains Transfac version information
    my $motif=Motif->new();      #motif
    while (<TRANSFAC>) {    ##Takes input from Transfac Matrix file(Parse file into array  Transfac_output)
            if (/^\/\//) {		## if line begins with // then this pushs the Matrix data onto array.
                                	## if the file begins with a // then it just ignores it and goes to the next line
        	    if (exists $motif->{POSITION_SCORING_MATRIX} ){
                        $motif->{VERSION}=$version;
                        $motif->{MOTIF_METHOD}="Transfac";
                        $motif->{FILE}=$filename;
			push @{$self->{MOTIFS}},$motif;
                        my $empty=Motif->new();
			$motif=$empty;		#Reset the variable for next matrix
                        }
		    else {next;}}
            elsif (!/^(XX|VV|P0|AC)/){         ##if it starts with AC, Delete AC and put the rest in Transfac_output[0], Should be Transfac table name  Example.  M00001
		s/(^..)\s+//;
                chomp $_;
                my $header=$transfac_header{$1};
                push (@{$motif->{$header}},$_);
            }
            elsif (/^AC/){  #If line begins with AC then parse to get identifier
                s/(^..)\s+//;
                chomp;
                $motif->{ACCESSION}=$_;
            }
            elsif (/^P0/) {##if it starts with PO, enter recursive loop to end of matrix
                my @transfac_matrix;
                my $consensus;
                while (<TRANSFAC>){
                    if (/^XX/){  ##if line begins with XX this signals the end of the matrix, exit loop with last command
                        last;}
                    else{       ##For all other lines in Matrix: split and put values corresponding to ACGT in matrix
                        my @line_matrix=split /\s+/,;
                        push @transfac_matrix, [@line_matrix[1..4]];
                        if (defined $line_matrix[5]){
                            $consensus.=$line_matrix[5]}        ## Append the letter values and add them to Transfac_output
                    }
                }
                $motif->{POSITION_SCORING_MATRIX}=\@transfac_matrix;
                $motif->{POSITION_PROBABILITY_MATRIX}=_create_ppm(\@transfac_matrix);
                $motif->{CONSENSUS}=$consensus;}
            elsif(/^VV/){
                s/(^..)\s+//;
                chomp $_;
                $version=$_;
            }
            else {next};
            }
    close TRANSFAC;
    return $self;
}


################################################################################
#
#   _create_ppm()
#
#   Description:  Creates Position_probability_matrix that is pseudo-counted
#                 a PWM can be calculated then by taking the log ratio against
#                 background frequency
#
#
sub _create_ppm{
    my ($psm,$background)=@_;
    my @psm=@{$psm};

    my @ppm = map {
            my @position=@{$_};
            my $total_counts;
            foreach (@position){
                $total_counts+=$_;
            }
            for(my $i=0;$i<4;$i++){
                $position[$i]=($position[$i]+(($total_counts)**(0.5))/4)/($total_counts+(($total_counts)**(0.5)));  #Pseudocount and calculate probability
                }
            \@position;} @psm;
    return \@ppm;  #return pseudocounted probability matrix
}


################################################################################
#
# _import_report
#
#   Description:  prints a report format of imported motifs
#
#   Work: make optional in import and new functions according to flag

sub _import_report {
    my ($self,$existing,$function,$filename)=@_;
    my @new_MotifS;
    foreach my $Motif (sort keys %$self){
        if (exists $existing->{$Motif}){
            next;
        }
        else {
            push @new_MotifS,$Motif;
        }
    }

    print "\n$function file: $filename\n";
    print "# of Imported Motifs: " . scalar(@new_MotifS) . "\n";
    printf "%-15s\t%-20s\t%-15s\n","ID\#", "ACCESSION", "SIZE OF PWM";
    print "--------------------------------------------------------\n";
    foreach my $Motif (@new_MotifS){
            my $ID=$self->{$Motif}->{FILE_ACCESSION}->[0];
            my $size=scalar(@{$self->{$Motif}->{POSITION_PROBABILITY_MATRIX}});
            printf "%-15s\t%-20s\t%-4d\n",$ID, $Motif, $size;
    }
}


##############################################################################################################
#
#   jaspar_import
#
#   Description:    Imports Jaspar SQL formated files into motifset using ANNOTATION.txt and DATA.txt files.
#                   The 3rd file INFO.txt is not used because they contain all the same data.
#                   If the values in the DATA section are integers then they are assigned to the
#                   POSITION_SCORING_MATRIX.  However, if numbers are decimal notation then they are
#                   put in POSITION_PROBABILITY_MATRIX.
#
#  Bug:  If previously input motifset, then numbering is wrong. when I go to add the PSM,PPM  Therefore, need to correct.
#
sub jaspar_import{
    my %default=("ANNOTATION"=>undef,
                 "DATA"=>undef
                 #"INFO"=>undef,
                 );
    my $self = shift;
    my %arg=(%default,@_);

    my %name_key;

    if ((!defined $arg{ANNOTATION}) || (!defined $arg{DATA})){
        die "Missing JASPAR filename";
    }

    if (exists $self->{Distance_Score}){
        delete $self->{Distance_Score};
    }

    open ANNOTATION, "< $arg{ANNOTATION}" || die "Couldn't open $arg{ANNOTATION}:$!\n";
    open PFM, "< $arg{DATA}" || die "Couldn't open $arg{DATA}:$!\n";
    #open INFO, "< $arg{INFO}" || die "Couldn't open $arg{INFO}:$!\n";

    my $motif;

    my $index= (exists $self->{MOTIFS}) ? scalar @{$self->{MOTIFS}} : 0 ;


    ##Import the annotations files from the ANNOTATION file
    while (<ANNOTATION>){
        chomp $_;
        if (/([A-Z]+\d+)\s(\w+)\s(.+)/){
            my $accession=$1;
            my $category=$2;
            my $info=$3;
            if (!exists $name_key{$accession}){
                $name_key{$accession}=$index;
                $index++;
                if (!defined $motif){
                    $motif=Motif->new();
                }
                else{
                    push @{$self->{MOTIFS}},$motif;
                        my $empty=Motif->new();
			$motif=$empty;		#Reset the variable for next matrix
                }
                $motif->{ACCESSION}=$accession;
                $motif->{MOTIF_METHOD}="JASPAR";
                push (@{$motif->{FILE}},($arg{ANNOTATION},$arg{DATA}));
            }

            if ($category eq "name"){
                $motif->{IDENTIFIER}=$info;}
            elsif ($category eq "medline"){
                $motif->{MEDLINE_ID}=$info;}
            elsif ($category eq "sysgroup"){
                $motif->{ORGANISM}.="\; $info";}
            elsif ($category eq "species"){
                $motif->{ORGANISM}=$info;}
            elsif ($category eq "consensus"){
                $motif->{CONSENSUS}=$info;}
            elsif ($category eq "class"){
                $motif->{PROTIEN_CLASS}=$info;
            }
            elsif ($category eq "acc"){
                my @acc;
                if ($info=~/\s,\s/g){
                    @acc=split(/,/,$info);
                    push (@{$motif->{PROTEIN}},@acc);
                }
                else {
                    $motif->{PROTEIN}=$info;
                }

            }
            elsif ($category eq "comment"){
                $motif->{COMMENTS}=$info;
            }
        }
    }
    push @{$self->{MOTIFS}},$motif;
    close ANNOTATION;

    ##Import the Position frequency matrices from the data file
    while (<PFM>){
        if (/([A-Z]+\d+)\s([ACGT])\s(\d+)\s([0-9.]+)/) { #integers
            my $accession=$1;
            my $nucleotide=$2;
            my $nucleotide_position= ($nucleotide eq "A") ? 0:
                                     ($nucleotide eq "C") ? 1:
                                     ($nucleotide eq "G") ? 2:
                                     3;
            my $position=$3;
            $position--;
            my $number=$4;
            my $matrix_type;
            my $motif_array=$name_key{$accession};
            if ($number=~/\d+\.\d+/){
                $matrix_type="POSITION_PROBABILITY_MATRIX";
            }
            else {
                $matrix_type="POSITION_SCORING_MATRIX";
            }

            $self->{MOTIFS}->[$motif_array]->{$matrix_type}->[$position]->[$nucleotide_position]=$number;
        }
    }

    close PFM;

    ##If necessary create POSITION_PROBABILITY_MATRIX for all imported motifs
    foreach my $array (values %name_key){
        if (exists $self->{MOTIFS}->[$array]->{POSITION_SCORING_MATRIX}){
            $self->{MOTIFS}->[$array]->{POSITION_PROBABILITY_MATRIX}=_create_ppm($self->{MOTIFS}->[$array]->{POSITION_SCORING_MATRIX});
        }
    }

    return $self;
}

###############################################################################################################
#
#   nestedmica_import
#
#   Description:    Imports Nested Mica XMS Motif file into Motifset.  Uses XMS identifiers such as <Motif>  to
#                  determine entries in Motif.  Nested Mica Motifs position probability matrix are pseudocounted by
#                  Nested Mica, therefore Motif's created by are not pseudocounted.  Nested Mica XMS lack counts
#                  therefore PSM are not computed. PWM are calculated using PPM and calculating log-likelihood
#                   Screens and issues warnings for negative values or columns that don't sum to 1.   Negative values
#                   changed to 0.0001 if the abs value sum is not 1.  (This accounts for mistakes where neg value used where in place of positive value.
#                   For columns that don't sum to one, the diffence is distributed among all columns.
#
#   Usage:  $Motifset->NestedMica_import("Nested Mica XMS filename")
#
#   Example:  $Motifset->NestedMica_import("Top-At-Motifs.xms")
#
#

sub nestedmica_import {
    my ($self,$filename)=@_;
    my $input_size=0;
    my $output_size=0;
    
    # # of motifs in self
    if (exists $self->{MOTIFS}){
        $input_size=$#{$self->{MOTIFS}} + 1;
    }
    
    #if Motifset has a Distance_Score matrix ; delete it.
    if (exists $self->{Distance_Score}){
        delete $self->{Distance_Score};
    }
    
    
    my $check_ppm=sub {
        my ($array,$name,$file,$column)=@_;
        my $sum=0;
        my $abs_sum=0;
        
        for (my $i=0;$i<4;$i++){
            $sum+=$array->[$i];
            $abs_sum+=abs($array->[$i]);
        }
        
        #if absolute value is close to one then change the sign
        if (abs(1-$abs_sum)<0.0001){
            for (my $i=0;$i<4;$i++){
                 if ($array->[$i]<0){
                    my $old_value=$array->[$i];
                    $array->[$i]=abs($array->[$i]);
                    print "Warning: Negative value of $old_value found in Motif:$name while processing File:$filename at column position: $column.  Value changed to $array->[$i]\n\n";
                 }
            }
            $sum = $abs_sum;
        }
        else {
            #check for negative values
            for(my $i=0;$i<4;$i++){
                if ($array->[$i]<0){
                    my $old_value=$array->[$i];
                    $array->[$i]=0.0001;
                    print "Warning: Negative value of $old_value found in Motif:$name while processing File:$filename at column position: $column.  Value changed to 0.0001\n\n";
                    #correct sum value
                    $sum=$sum-$old_value+0.0001;
                }    
            }
        }
        
        
        if (abs(1-$sum)>0.0001){
            my $difference=1-$sum;
            my @old_array=@$array;
            if ($difference<-0.0001){
                for(my $i=0;$i<4;$i++){
                    $array->[$i]-=$difference/4;
                }
            }
            elsif ($difference>0.0001){
                for(my $i=0;$i<4;$i++){
                    $array->[$i]+=$difference/4;
                }
            }
            print "Warning: Column sum error: Sum of $old_array[0] + $old_array[1] + $old_array[2] + $old_array[3] = $sum\n\n";
            print "Column $column of Motif: $name in File:$filename adjusted to equal 1 by a value of $difference \n\n";
        }
        return $array;
    };
    
    
    # Control variables to make sure imports were done correctly
    my @headers=(0,0,0,0,0,0,0,0,0,0,0);
    #0my $top_motifset_header=0;
    #1my $bottom_motifset_header=0;
    #2my $top_motif_header=0;
    #3my $bottom_motif_header=0;
    #4my $accession_header=0;
    #5my $top_weightmatrix_header=0;
    #6my $bottom_weightmatrix_header=0;
    #7my $top_column_header=0;
    #8my $bottom_column_header=0;
    #9my $weights_header=0;
    #10my $threshold_header=0;
    
    my $motif=Motif->new();
    open (NESTED, "<$filename") || die "Couldn't Open $filename:$!\n;";
    while (<NESTED>){
        if (/\<motifset/){
            $headers[0]++;
        }
        elsif (/\<motif\>/){
            $headers[2]++;}
        elsif (/\<name\>(.+)\<\/name\>/){
            $headers[4]++;
            #print $1;
            $motif->{ACCESSION}=$1;}
        elsif (/\<weightmatrix/) {
            $headers[5]++;
            my @ppm;
            my @ppm_position;
            while (<NESTED>){
                if (/\<column/){
                    $headers[7]++;
                }
                elsif (/\<\/weightmatrix\>/){
                    $headers[6]++;
                    last;}
                elsif (/\<\/column\>/){
                    $headers[8]++;
                    &$check_ppm(\@ppm_position,$motif->{ACCESSION},$filename,scalar @ppm);
                    push(@ppm, [@ppm_position]);
                    undef @ppm_position;}
                elsif (/\<weight symbol="([a-z]+)"\>(\d\.\d+|\d|-\d\.\d+|-\d)/){
                    $headers[9]++;
                    if      ($1 eq "cytosine")  {$ppm_position[1]=$2;}
                    elsif   ($1 eq "thymine")   {$ppm_position[3]=$2;}
                    elsif   ($1 eq "guanine")   {$ppm_position[2]=$2;}
                    elsif   ($1 eq "adenine")   {$ppm_position[0]=$2;}
                     }
                }
            $motif->{POSITION_PROBABILITY_MATRIX}=\@ppm;
        }
        elsif (/\<threshold\>(\d\.\d+|\d)\<\/threshold\>/){
            $headers[10]++;
            $motif->{THRESHOLD}=$1;
        }
        elsif (/\<\/motif\>/){
            $headers[3]++;
            $motif->{FILE}=$filename;
            push @{$self->{MOTIFS}},$motif;
            my $empty=Motif->new();
            $motif=$empty;
        }
        elsif (/\/motifset/){
            $headers[1]++;
        }
        
    }
    
    
    #check for import errors
    my $import_errors=undef;
    if ($headers[0]!=$headers[1])   { $import_errors.="Problem with motifset tags\n";}
    if ($headers[2]!=$headers[3]){$import_errors.="Problem with motif tags\n";}
    if ($headers[4]!=$headers[2]){$import_errors.="# of name tags doesn't match # of motif tags\n";}
    if ($headers[5]!=$headers[6]){$import_errors.="Problem with weightmatrix tags\n";}
    if ($headers[7]!=$headers[8]){$import_errors.="Problem with column tags\n";}
    if ($headers[9]/4!=$headers[7]){$import_errors.="Problem with weight tags\n";}
    if ($headers[10]!=$headers[5]){$import_errors.="Problem with threshold tags\n";}
    
    if (defined $import_errors){
        print "Formatting errors found for $filename, not all motifs imported or imported improperly\n\n";
        print $import_errors;
    }
    
    if (exists $self->{MOTIFS}){
         $output_size=$#{$self->{MOTIFS}} + 1;
    }
    if ($input_size==$output_size){
        print "Error: No motifs were imported from $filename\n";
    }
    else {print $output_size-$input_size . " motifs imported from $filename.   Number of motifs in motifset is $output_size\n";}
    
    
    close NESTED;
    return $self;
}


###################################################################################################
#  meme_Import
#
#  Description:    Meme_Import imports MEME Motifs text file into Motifset. Motifset created by parsing
#                  line headers of the text report.  PSM calculated using Position Probability Matrix
#                  and Motif_Locations at the end of the Motif creation.   PWM are then calculated by
#                  pseudocounting the PSM and calculating log likelihood ratio.
#
#  Note:  This will not work to parse MEME HTML report.
#
#  Usage:  $motifset->meme_import("Meme Text report file")
#
#  Example:  $motifset->meme_import("Meme Output DNA.txt")
#
#
sub meme_import {
    my ($self,$filename)=@_;
    open (MEME, "< $filename") || die "Couldn't Open $filename:$!\n;";

    if (exists $self->{Distance_Score}){
        delete $self->{Distance_Score};
    }

    my $version;
    my $datafile;
    my $alphabet;
    my $job;
    my @sequences;
    my $commandline;
    my $background_freq=SEQ_FREQ->new();
    my $dataset_freq=SEQ_FREQ->new();
    while (<MEME>){
        if (/MEME job (\d+)/){
            $job=$1;}
        elsif (/MEME (version [0-9.]+)/){
            $version=$1;}
        elsif (/DATAFILE=(.+)/){
            $datafile=$1;}
        elsif (/ALPHABET=(.+)/){
            $alphabet=$1;}
        elsif (/Sequence name/){
            while(<MEME>){
                if (/\-+/){
                    next;}
                elsif (/\*+/){
                    last;}
                else {
                    my @line=split(/\s+/,$_);
                    push @sequences,$line[0];
                    if (defined $line[3]){
                        push @sequences,$line[3];}
                }
            }
        }
        elsif (/command: (.+)/){
            $commandline=$1;}

        #import sequence frequency for dataset
        elsif (/Letter frequencies in dataset:/){
            $_=<MEME>;
            if (/A [0-9.]+ /){
                my @freq_line=split(/\s+/,$_);
                $dataset_freq->{1}={"A" => $freq_line[1],
                            "C" => $freq_line[3],
                            "G" => $freq_line[5],
                            "T" => $freq_line[7]};
            }
            else {
                    die "meme_import: letter frequencies not standard nucleotide";}
        }

        #import background letter frequencies for
        elsif (/Background letter frequencies/){
            $_=<MEME>;
            if (/A [0-9.]+ /){
                my @freq_line=split(/\s+/,$_);
                $background_freq->{1}={"A" => $freq_line[1],
                            "C" => $freq_line[3],
                            "G" => $freq_line[5],
                            "T" => $freq_line[7]};
            }
            else {
                    die "meme_import: background frequencies not standard nucleotide";}
        }
        elsif (/\s+Motif (\d+) Description/){
                my $motif=Motif->new();
                $motif->{ACCESSION} = "Motif $1";
                $motif->{VERSION}= $version;
                $motif->{DATAFILE}=$datafile;
                $motif->{COMMAND_LINE}=$commandline;
                $motif->{ALPHABET}= $alphabet;
                $motif->{JOB}= $job;
                $motif->{SEQUENCES}=\@sequences;
                $motif->{BACKGROUND_FREQUENCY}=$background_freq;
                $motif->{DATASET_FREQUENCY}=$dataset_freq;
                $motif->{FILENAME}=$filename;
                while (<MEME>){
                    if (/Time/){
                        last;}
                    elsif (/\s+Motif \d+ sites sorted by position p-value/){
                        $_=<MEME>;
                        while (<MEME>){
                            if (/Sequence name/ || /^\-{13}\s+/){
                                next;}
                            elsif (/(\S+)\s+([+|-| ])\s+(\d+)\s+[0-9.e-]+/){
                                my @motif_pos=($1,$2,$3);
                                push (@{$motif->{MOTIF_LOCATIONS}},\@motif_pos)
                                #$motif->{MOTIF_LOCATIONS}set_attribute("ML",\@Motif_pos);
                                }
                            elsif (/^\-{80}/){
                                last;}
                            else {
                                next;}}}

                    #import position probability matrix
                    elsif (/\s+Motif \d+ position\-specific probability matrix/){
                        $_=<MEME>;
                        my @ppm;
                        while(<MEME>){
                            if (/(\s+[0-9.]+){4}/){
                                my @row=split /\s+/, $_;
                                push(@ppm,[@row[1..4]]);}
                            elsif (/^\-{80}/){
                                last;}
                            else {
                                next;}
                            }

                        #convert Position probability matrix to Position scoring matrix
                        my @psm;
                        my $total_counts=scalar (@ppm);
                        @psm = map { my @position=map $_*$total_counts,@{$_};
                                    \@position;} @ppm;
                        $motif->{POSITION_SCORING_MATRIX}=\@psm;

                        #create pseudo_counted position probability matrix
                        $motif->{POSITION_PROBABILITY_MATRIX}=_create_ppm(\@psm);
                        }

                    #import regular expressions for motif
                    elsif (/\s+Motif \d+ regular expression/){
                        my $regex;
                        while (<MEME>){
                            if (/\-+/){
                                next;}
                            elsif (/[ACGT\[\]]+/){
                                $regex=$_;}
                            else{
                                last;}
                            }
                        chomp $regex;
                        $motif->{REGULAR_EXPRESSION}=$regex;}
                    else {
                        next;}
                    }

                #import into Motifset;
                push @{$self->{MOTIFS}},$motif;
            }
        else {
            next;}
    }
    close MEME;
    return $self;
}



##################################################################################################
#  weeder_import
#
#  Description:    Import Weeder text output .wee format in Motifset.   Uses line headers to parse the weeder
#                  file into Motif.   PPM and PWM are calculated from the All Occurences Frequency Matrix.
#                  PWM calcuted from All Occurences. See PWM_matrix.
#
#  Usage:  $Motifset->Weeder_import("Weeder wee text file")
#
#  Example:  $Motifset->Weeder_import("G1S.cycle.fasta.wee")
#
#
sub weeder_import {
    my ($self,$filename)=@_;
    open (WEEDER, "< $filename") || die "Couldn't Open $filename:$!\n;";

    #delete DISTANCE_SCORE
    if (exists $self->{Distance_Score}){
        delete $self->{Distance_Score};
    }

    my @commandlines;
    my @top_motifs;
    my @sequences;
    my @motifs;
    my $organism;
    my $i=0;

    while (<WEEDER>){
        if (/^Organism code: (\w+)/){
            $organism =$1;}

        #store commandline in array
        elsif (/Searching for motifs/){
            my @command_info;
            while (<WEEDER>) {
                if (/^\.\/(.*)\)/){
                    $command_info[0]=$1;}
                elsif (/\d{1,2}\).([ACGT]{5,15})/){
                    push(@{$command_info[1]},$1);}
                elsif (/Elapsed time/){
                    push @commandlines,\@command_info;
                    last;}}}
        #store sequence accessions in sequence array
        elsif (/^Sequence \d+\ : \>(\S+)/){
            push(@sequences,$1);
        }

        #starting of individual motif
        elsif (/\*+ MY ADVICE \*+/){
            my $motif=Motif->new();
            while(<WEEDER>){

                #Motif and redundant motifs section
                if (/(^[ACGT]{5,15})\n/){
                    push @{$motif->{MOTIF_SEQUENCES}},$1;
                    $motif->{ACCESSION}=$1;
                    foreach (@commandlines){
                        if (grep (/$1/,@{$$_[1]})){
                            $motif->{COMMAND_LINE}=$$_[0];
                            last;}
                        else {next}
                    }
                }
                elsif (/[ACGT]{5,15}.-/){
                    chomp;
                    s/\s*//;
                    my @motifs=split /-/,$_;
                    map {push @{$motif->{MOTIF_SEQUENCES}},$_} @motifs;
                }

                #best occurrences motifs
                elsif (/^Best occurrences\:/){
                    my @motif_positions;
                    while(<WEEDER>){
                        if (/\d+ [+|-]/){
                            my @temp_row=split(/\s+/);
                            my $sequence=$temp_row[1];
                            $sequence--;
                            $temp_row[1]=$sequences[$sequence];
                            @temp_row=@temp_row[1,2,3,4,5];
                            push @{$motif->{MOTIF_LOCATIONS}}, \@temp_row;}
                        elsif (/^\s+/){
                            last;}}}
                elsif (/\s+Frequency Matrix/){
                    my @all_occurrences;
                    my @best_occurrences;
                    while (<WEEDER>){
                        if (/^\d+/){
                            my @row=split(/\s+/);
                            push(@all_occurrences,[@row[1,2,3,4]]);
                            push(@best_occurrences,[@row[5,6,7,8]]);}
                        elsif (!/[A-Za-z0-9]/){
                            last;
                            }
                        }
                    $motif->{POSITION_SCORING_MATRIX}=\@all_occurrences;
                    $motif->{BEST_OCCURENCES_SCORING_MATRIX}=\@best_occurrences;
                    $motif->{POSITION_PROBABILITY_MATRIX}=_create_ppm(\@all_occurrences)
                }

                elsif (/\=+/){
                    $motif->{FILENAME}=$filename;
                    $motif->{ORGANISM}=$organism;
                    push @{$self->{MOTIFS}},$motif;
                    $motif=Motif->new();
                }
                else {
                    next;
                }
            }
        }
    }
    close WEEDER;
    return $self;
}

##################################################################################################
#
#   transfac_output
#
#   Description:    Outputs the Motifset to a file in a transfac compatible format .   However, will not work
#                  on a Motifset created from Nested Mica XMS file, because Nested Mica does not
#                  contain a Position Scoring Matrix.
#
#  Note:   Need to work out a fix for Nested Mica XMS file using Seq_Score  or User input of # of
#          Motifs found.  Thereby, allowing for the creation of a scoring matrix
#
#  Usage:  $Motifset->transfac_output("My output file.txt")
#
#   Work: Standarize outputs and possiblity of outputing transfac-like output instead
#           of only transfac compatibility
#


sub transfac_export{
    my ($self,$filename)=@_;
    my @Categories=("ACCESSION",
                    "IDENTIFIER",
                    "DATE_AUTHOR_UPDATE_DATES",
                    "COPYRIGHT",
                    "NAME_OF_BINDING_SITE",
                    "SHORT_FACTOR_DESCRIPTION",
                    "LIST_OF_LINKED_FACTOR_ENTRIES",
                    "POSITION_SCORING_MATRIX",
                    "STATISTICAL_BASIS",
                    "BINDING_SITES_UNDERLYING_MATRIX",
                    "COMMENTS",
                    "MEDLINE_ID",
                    "REFERENCE_#",
                    "REFERENCE_AUTHORS",
                    "REFERENCE_TITLE",
                    "REFERENCE_DATA");
    my %prompts =("ACCESSION" => "AC",
                  "FILE_ACCESSION" => "FA",
                  "IDENTIFIER" => "ID",
                  "DATE_AUTHOR_UPDATE_DATES"=> "DT",
                  "COPYRIGHT"=>"CO",
                  "NAME_OF_BINDING_SITE" => "NA",
                  "SHORT_FACTOR_DESCRIPTION" => "DE",
                  "LIST_OF_LINKED_FACTOR_ENTRIES" => "BF",
                  "POSITION_SCORING_MATRIX"=>"PO",
                  "STATISTICAL_BASIS"=>"BA",
                  "BINDING_SITES_UNDERLYING_MATRIX"=>"BS",
                  "COMMENTS"=>"CC",
                  "MEDLINE_ID"=>"RX",
                  "REFERENCE_#"=>"RN",
                  "REFERENCE_AUTHORS"=>"RA",
                  "REFERENCE_TITLE"=>"RT",
                  "REFERENCE_DATA"=>"RL");
    if (defined $filename){
        open OUT, "> $filename";
        select (OUT);}
    foreach my $Motif (keys %$self){
        if (ref $self->{$Motif} eq "Motif"){
            foreach my $key (@Categories){
                if (exists $self->{$Motif}->{$key}){
                    if (ref $self->{$Motif}->{$key} eq "array"){
                        foreach my $entry (@{$self->{$Motif}->{$key}}){
                            print "$prompts{$key}\t$entry\n";
                        }
                    }
                    elsif (ref $self->{$Motif}->{$key} eq "PSM_MATRIX"){
                        print "PO\tA\tC\tG\tT\n";
                        my $n="00";
                        foreach my $counts (@{$self->{$Motif}->{$key}}){
                            my @counts=@$counts;
                            printf "%02s\t%d\t%d\t%d\t%d\n",$n,$counts[0],$counts[1],$counts[2],$counts[3];
                            $n++;
                        }
                    }
                    elsif ((ref $self->{$Motif}->{$key} eq "PPM_MATRIX")&&(!exists $self->{$Motif}->{POSITION_SCORING_MATRIX})){
                        my $output=select STDOUT;
                        print "Position Scoring Matrix doesn't exist for Motif.   Need to generate PSM_Matrix using sequence_set\n";
                    }
                    else {
                    print "$prompts{$key}\t$self->{$Motif}->{$key}\n";
                    }
                    print "XX\n";
                }
            }
            print "\/\/\n";
        }
        else {next;}
    }
    $|=1;
    close OUT;
    select (STDOUT);
}


#################################################################################################
#
#   nestedmica_output
#
#   Description:  output motifset to XMS file format
#
sub nestedmica_export{
    my ($self,$filename)=@_;
    my @template=@{main::nm_template()};
    
    open (XMS, ">$filename.xms") || die "Couldn't open file \n";
    
    print XMS "$template[0]\n"; #print header for xms file
    foreach my $motif (@{$self->{MOTIFS}}){
        my $name=$template[6];
        $name=~s/\$motif/$motif->{ACCESSION}/g;
        
        my $columns=scalar @{$motif->{POSITION_PROBABILITY_MATRIX}};
        my $columns_output=$template[7];
        $columns_output=~s/\$columns/$columns/g;
        
        print XMS "$template[5]\n";  #print start of motif
        print XMS "$name\n";
        print XMS "$columns_output\n";
        my $i=0;
        foreach my $column (@{$motif->{POSITION_PROBABILITY_MATRIX}}){
            my $column_header=$template[8];
            my $a=$template[9];
            my $c=$template[10];
            my $g=$template[11];
            my $t=$template[12];
            $column_header=~s/\$column/$i/g;
            $a=~s/\$a_value/$column->[0]/g;
            $c=~s/\$c_value/$column->[1]/g;
            $g=~s/\$g_value/$column->[2]/g;
            $t=~s/\$t_value/$column->[3]/g;
            print XMS "$column_header\n$a\n$c\n$g\n$t\n$template[13]\n";
            $i++;
        }
        
        my $threshold=$template[15];
        if (defined $motif->{THRESHOLD}){
            $threshold=~s/\$threshold/$motif->{THRESHOLD}/g;
        }
        else{
            $threshold=~s/\$threshold/0/g;
        }
        
        print XMS "$template[14]\n$threshold\n$template[16]\n";
    }
    print XMS "$template[17]\n";
    
    $|=1;
    close XMS;

}




#################################################################################################
#
#   randomize_motifset
#
#   Description:  Randomize motif based on column positions (Kluth Shuffle)
#                 Returns a motifset.
#   Usage:  $new_motifset->randomize_motifset("MOTIFSET"=>$motifset,"MOTIF_SELECTION"=>[1,3,4,9],"VERBOSE"=>1);
#           This will take $motifset and extract motifs found in the motifset's array 1,3,4,9 and rearrange columns randomly.
#           Output will be a motifset with 4 motifs that are randomized.   If verbose is defined as anything, then output to
#           console the motif name and new arrangement of columns.
#
#   Fix:  CHANGE PRINT VERBOSE TO PRINTF


sub randomize_motifset {
    #set default values
    my %default=("MOTIFSET"=>undef,
                 "MOTIF_SELECTION"=>undef,
                 "VERBOSE"=>undef);
    my $new=shift;
    my %arg=(%default,@_);
    
    #check to make sure that MOTIFSET is defined and is a Motifset object
    if (!defined $arg{MOTIFSET} || ref $arg{MOTIFSET} ne "Motifset"){
        die "MOTIFSET was not defined or not a Motifset object\n";
    }
    
    #Check to make sure that the MOTIF_SELECTION is an array if it is defined
    if (defined $arg{MOTIF_SELECTION} && ref $arg{MOTIF_SELECTION} ne "ARRAY" ){
        die "MOTIF_SELECTION must be a reference to an array"
    }
    
    my $old=$arg{MOTIFSET};
    
    if (defined $new->{MOTIFS}){
        delete $new->{MOTIFS}
    }
    
    if (defined $arg{VERBOSE}){
        printf "Motif\t\t\t\tNew Order\n";
        print "--------------------------------------------------\n"
    }
    
    if ((defined $arg{MOTIF_SELECTION}) && (ref $arg{MOTIF_SELECTION} eq "ARRAY")){
        foreach my $motif (@{$arg{MOTIF_SELECTION}}){
            
            #Check to see if $motif is valid selection from Motifset
            if (!exists $old->{MOTIFS}->[$motif]){
                my $size=scalar @{$old->{MOTIFS}};
                die "$motif is out of range.  There are only $size motifs in the motifset\n";
            }
            
            my $accession=$old->{MOTIFS}->[$motif]->{ACCESSION};
            if (defined $arg{VERBOSE}){
                printf "$accession\t\t";
            }
            my @pwm=main::randomize_pwm($old->{MOTIFS}->[$motif]->{POSITION_PROBABILITY_MATRIX},$arg{VERBOSE});
            my $new_motif=Motif->new();
            $new_motif->{ACCESSION}="Randomized " . $accession;
            $new_motif->{POSITION_PROBABILITY_MATRIX}=\@pwm;
            push @{$new->{MOTIFS}},$new_motif;
        }
    }
    else{       #if selection not given, assume that whole motifset is to be randomized.
        foreach my $motif (@{$old->{MOTIFS}}){
            my $accession=$motif->{ACCESSION};
            if (defined $arg{VERBOSE}){
                print "$accession\t\t";
            }
            my @pwm=main::randomize_pwm($motif->{POSITION_PROBABILITY_MATRIX}, $arg{VERBOSE});
            my $new_motif=Motif->new();
            $new_motif->{ACCESSION}="Randomized " . $accession;
            $new_motif->{POSITION_PROBABILITY_MATRIX}=\@pwm;
            push @{$new->{MOTIFS}},$new_motif;
        }
    }
    
    return $new;
}



#----------------------------------------------------------------------------------------------------------------------------------------------------

##################################################################################################
#  SEQ_FREQUENCY PACKAGE
#
#  Description:    SEQ_FREQUENCY package used by SEQUENCE and SEQUENCE_SET function to calculate
#                  freqeuncy table for wordsizes.  SEQ_FREQUENCY is a hash, where HASH OF HASH.
#                  Stored in format of {Wordsize}->{Word}->Count.  Example  {'3'}->{AAG} => 3
#
#   

package SEQ_FREQ;

##new - accessed from SEQUENCE or SEQUENCE_SET functions
sub new{
    my $class=shift;
    my $self=bless {},$class;
    return $self;
}

##################################################################################################
#  frequency
#
#  Description:  Counts frequency of word occurences within SEQUENCE or SEQUENCE_SET object.
#                Default options are to compute the frequency for wordsize 1,2,3 using forward and
#                reverse sequence.   However, if the user wants to only calculate the forward or reverse strand for
#                a different size can supply separate flag.
#
#
#   Default 1,2,3 nucleotide frequencies
#   Flag
#   -n, where n=# computes n-word size frequency
#   -f,-r,-b forward sequence only, reverse seqence only, both forward and reverse (Default: Forward);
#
#
sub frequency{
    my ($self,$seq,$n_word,$seq_flag)=@_;
    my $sequence=" ";
    my $size;
    print ref($seq);

    if (!defined $n_word){
        $n_word=3;
    }

    #If the only a sequence is handed to function then depending on the flags
    #it will use or append the reverse complement to the sequence.
    if ((ref($seq) eq "SCALAR") || !(ref($seq))){
        my $rev_comp=$seq;
        $rev_comp=~tr/ACGTN/TGCAN/;
        $rev_comp=reverse($rev_comp);

        if (defined $seq_flag && $seq_flag eq "\-r"){
            $sequence=$rev_comp;
            $size=length $sequence;}

        elsif (defined $seq_flag && $seq_flag eq "\-b"){
            $sequence=" " . $seq . " " . $rev_comp;
            $size=length $sequence;}

        else{
            $sequence=$seq;
            $size=length $sequence;}
    }

    #If an array of sequences are handed to the function.  It will append the
    #sequences together with intervening spaces, so the sequences don't generate words
    #that don't exists within the sequences
    elsif (ref $seq eq "ARRAY"){
        foreach my $seqs (@$seq){
            my $rev_comp=$seqs;
            $rev_comp=~tr/ACGTN/TGCAN/;
            $rev_comp=reverse($rev_comp);
            if (defined $seq_flag && $seq_flag eq "\-r"){
                $sequence.=" " . $rev_comp;}
            elsif (defined $seq_flag && $seq_flag eq "\-b"){
                $sequence.=" " . $seqs . " " . $rev_comp;}
            else{
                $sequence.= " " . $seqs;}
        }
        $size=length $sequence;
    }

    #compute the frequency table for all word sizes less than n
    for (my $n=1;$n<=$n_word;$n++){
        my %counts;
        my $skipped=0;

        #counts the words unless they have a space, in which case it increments
        # skipped and goes to next word.   This is so it won't count words that
        # aren't in the sequences when they were appended to each other.
        for (my $i=0;$i<$size-($n-1);$i++){
            my $substr=substr($sequence,$i,$n);
            if ($substr=~/\s+/){
                $skipped++;
                next;}
            else {

                $counts{$substr}++;}
        }

        #Calculated frequency from counts
        my $divisor=$size-($n-1)-$skipped;
        foreach my $key (keys %counts){
            $counts{$key}/=$divisor;
         }
        $self->{$n}=\%counts;
    }
    return $self;
}



##########################################################################################################################################

package main;

##_query_user queries the user as to which of the sequences or Motifs in a set to use for function that takes only Motif or sequence
sub _query_user {
    my $self=shift;
    my $output;
    my @query;
    my @index;
    if (ref $self eq "Motifset"){ #Determine if Motifset or Sequence_set and then set @query to the appropriate queues
        @query=("Motif", "Motif");
    }
    else { print "_query_user doesn't handle" . ref $self;
          return 1;
          }

    SELECTION:
    my $n=1;
    foreach my $keys (sort %$self){  #print and push all availabe Motifs or sequences to @index
        if (ref $self->{$keys} eq $query[0]){
            print "$n) $keys\n";
            push @index, $keys;
            $n++;
        }
    }
    print "Which $query[1] would you like (Example: 1) :" ; #Prompt user to select one of the previously displayed Motifs or sequences
    my $selection=<STDIN>;
    chomp $selection;
    if (!($selection=~/\d+/) || $selection>=$n){   #Warn if user input doesn't correspond to listed options :  Repeat query of user
        print "Your selection doesn't exist: reselect a $query[1]\n";
        goto SELECTION;
    }
    else{   #Return the user selected Motif or sequence
        $selection-=1;
        $output=$self->{$index[$selection]};
        return $output;
    }
}

##simple_generate ($sequence,size);
#Generate sequence using mono-nucleotide frequency
#
#
#Usage:  $sequence->generate_sequence($old_sequence,sequence_size);
#
#Example:  $sequence->generate_sequence($old_sequence,1000);
#
#
sub simple_generate{
    my ($seq_freq,$size)=@_;
    my $template;

    my $frequency=$seq_freq->{'1'};
    my $seq;
    foreach my $key (keys %$frequency){
        my $freq=$frequency->{$key};
        my $n=$size*$freq;
        use integer;
        $n/=1;
        for (my $i=0;$i<$n;$i++){
            $seq.=$key;
        }
    }
    my $n=0;
    my @array=sort {$frequency->{$b} <=> $frequency->{$a}} keys %$frequency;
    while (length $seq < $size){
        if ($n==4){
            $n=0;}

        $seq.=$array[$n];
        $n++;
    }

    my @seq=split('',$seq);

    #Randomize using kluth shuffle
    for ( my $i = $#seq; $i > 0; --$i ) {
        my $j = int rand( $i + 1 );
        @seq[ $i, $j ] = @seq[ $j, $i ];
    }
    $seq=join('',@seq);

    return $seq;
}

##
# generate_sequence
# Differs from quick_generate_sequence by using an table to make the distribution more accurate
#   For example: if TTT is used then it is less likely to be used in the next round of generation and the remaining
#               words TTA, TTG, TTC would have increased probabilities in the next run
#   First $wordsize-1 word generated without table;
#   Therefore, use simple generate for more accurate sequence length
#   Also, if word is shorter than orginal sequence then likely that words are missing from the new sequence.
#   Returned sequence is only approximately the asked for length, as the table may not be able to figure out the best solution to reach the exact length.
#
# Works fine for up to wordsize 7  but if asked for 10,000 and wordsize 7--returns sequence of approximately 6000
#
# USAGE STATEMENT:
#   usage:  $sequence->acc_sequence($old_sequence,sequence_length,wordsize)
#   Example:  $sequence->acc_sequence($old_sequence,9000,5);
#
#               $old_sequence: can be SEQUENCE object or SEQUENCE_SET object
#
#  As stated above, the generated sequence may be significantly less than the asked for length.   If accurates size is needed, use the simple quick_sequence_generate.
#
sub generate_sequence{
    my ($seq_freq,$length,$wordsize)=@_;
    my $sequence;

    if ($wordsize==1){
        $sequence=simple_generate($seq_freq,$length);
        return $sequence}


    $length=$length+200;
    my $max_length=0;

    if (!exists $sequence->{SEQ_FREQ}->{$wordsize}){
        die "Seq_Freq doesn't have necessary wordsize $wordsize\n";
    }

    my %seq_freq=%{$seq_freq};
    my %table;

    foreach my $word (sort keys %{$seq_freq{$wordsize}}){
        my $value=sprintf("%d",$seq_freq{$wordsize}{$word}*$length);
        if ($wordsize>1 && $value>0){
            $table{$wordsize}{$word}=$value;
            my $letter=substr($word,0,$wordsize-1);
            $table{$wordsize-1}{$letter}+=$table{$wordsize}{$word}
        }
    }
    my $new_seq="";
    for (my $i=0;$i<$length;$i++){
        my $random=rand(10);
        if ($i==0 || $wordsize==1){
            my $cumulative_total=0;
            foreach my $letter (sort keys %{$seq_freq{1}}){
                $cumulative_total+=10*$seq_freq{1}{$letter};
                if ($random<=$cumulative_total){
                    $new_seq.=$letter;
                    last;
                }
            }
        }
        elsif ($i>0 && $i<$wordsize-1){
            my $last=substr($new_seq,0,$i);
            my $cumulative_total=0;
            foreach my $word (sort keys %{$seq_freq{$i+1}}){
                if ($word=~/^$last/){
                    $cumulative_total+=10*($seq_freq{$i+1}{$word})/($seq_freq{$i}{$last});
                    if ($random<=$cumulative_total){
                        $new_seq.=substr($word,$i,1);
                        last;
                    }
                }
                else {next;}
            }

        }
        elsif ($i>=$wordsize-1){
            my $new_length=length $new_seq;
            my $index=$new_length-($wordsize-1);
            my $last=substr($new_seq,$index,$wordsize-1);
            #if ($last eq "AA"){
            #    print "here we go again";}  Only used to debug problem .   Can delete when finished
            my $cumulative_total=0;
            EVALUATE:
            foreach my $word (sort keys %{$table{$wordsize}}){
                if ($word=~/^$last/){

                    if ($table{$wordsize-1}{$last}==0){
                        my $seq_length=length $new_seq;
                        if ($max_length<$seq_length){
                            $max_length=$seq_length}
                        if ($max_length-$seq_length>=10){
                            $i=$length;
                            last;
                        }
                        $new_seq=substr($new_seq,0,$seq_length-1);
                        $seq_length=length $new_seq;
                        last;
                    }


                    $cumulative_total+=10*($table{$wordsize}{$word}/$table{$wordsize-1}{$last});
                    #print "Last : $last   Word: $word  CT:  $cumulative_total   Random:  $random\n";  #used for debug can delete when finished

                    if ($random<=$cumulative_total || abs($random-$cumulative_total)<0.1){  #difference b/c doesn't always add upto 10 and would produce error
                        my $temp_seq= $new_seq . substr($word,$wordsize-1,1);
                        my $temp_len=length $temp_seq;
                        my $index=$temp_len-($wordsize-1);
                        my $next_last=substr($temp_seq,$index,$wordsize-1);
                        if (!exists $table{$wordsize-1}{$next_last}){
                            $random=rand(10);
                            $cumulative_total=0;
                            goto EVALUATE;
                        }
                        else{
                            $new_seq.=substr($word,$wordsize-1,1);
                            $table{$wordsize}{$word}--;
                            $table{$wordsize-1}{$last}--;
                            last;
                        }
                    }
                }
                else {next;}
            }
        }
    }

    return $new_seq;
}



##############################################################################
#   randomize($sequence)
#
# Uses Knuth shuffle algorithm to randomize the sequence
# Code used from #http://www.perlmonks.org/?node_id=441060
# Returns new shuffled sequence object
#
# Usage $seq=randomize($sequence)
#

sub randomize {
    my ($sequence)=@_;
    if (!defined $sequence){
        die "No SEQUENCE given\n";
    }

    my @seq=split('',$sequence);
    for ( my $i = $#seq; $i > 0; --$i ) {
        my $j = int rand( $i + 1 );
        @seq[ $i, $j ] = @seq[ $j, $i ];
    }
    $sequence=join('',@seq);

    return $sequence;
}


################################################################################
#   simple_shuffle($sequence,n_shuffes,wordsize)
#   Performs a riffle (shuffle like deck of cards) on sequence
#
#   Usage: $newseq=simple_shuffle($sequence,n_shuffles,n_wordsize)
#
#   Example: $newseq=simple_shuffle($sequence,4,1);
#       Shuffles the sequence 4 times using single nucleotide as basic unit
#
#   Can use Sequence_set in place of $sequence.    Will query user for sequence.
#
#   If shuffle number or word size not chosen, then number 1 chosen for each parameter.

sub simple_shuffle {
    my ($sequence,$n_shuffles,$n_wordsize)=@_;

    if (!defined $sequence){
        die "No SEQUENCE given to function.";
    }

    if (!defined $n_shuffles) {
        $n_shuffles=1;
    }
    if (!defined $n_wordsize){
        $n_wordsize=1;
    }

    my $size=length $sequence;

    #Shuffles $n_shuffles amount of times
    for (my $n=0;$n<$n_shuffles;$n++){
        #Split sequence into halfs
        my $first_half=substr($sequence,0,$size/2);
        my $second_half=substr($sequence,$size/2,$size);
        my $new_seq="";
        #Concatenate $n_wordsize substring to new string; first from second half and then from second half
        for (my $i=0;$i<length $second_half;$i+=$n_wordsize){
            my $substr1=substr($second_half,$i,$n_wordsize);
            my $substr2=substr($first_half,$i,$n_wordsize);
            $new_seq.=$substr1 . $substr2;}
        #$new_seq then re-assigned to $seq for further shuffling
        $sequence=$new_seq;
    }

    return $sequence;
}

###############################################################################
#   generate_motif
#
#   Description:   Generates nucleotide motif sequence based on motif object and background frequency
#
#
sub generate_motif {
    my ($motif,$background)=@_;
    my %back_freq;

    if (defined $background){
        %back_freq=%{$background->{1}};
    }


    if (ref $motif eq "Motifset"){
        $motif=main::_query_user($motif);}

    my @motif=@{$motif->{POSITION_PROBABILITY_MATRIX}};
    my $motif_name=$motif->{ACCESSION};
    my $motif_size=scalar @motif;
    my $motif_seq="";
    foreach my $positions (@motif){
        my $rand=rand(10);
        my $A=10*@$positions[0];
        my $C=(10*@$positions[1])+$A;
        my $G=(10*@$positions[2])+$C;
        $motif_seq.= ($rand<$A) ? "A" :
                     ($rand<($C)) ? "C" :
                     ($rand<($G)) ? "G" :
                     "T";
        }

    my $LL; #log likelihood
    for (my $i=0;$i<length $motif_seq;$i++){
        my $letter=substr($motif_seq,$i,1);
        $LL+=   log($motif[$i][0]/$back_freq{$letter})/log(2);
    }

    return $motif_seq;
}



#################################################################################################
#
#   randomize_pwm
#
#   Description:  Randomize motif based on column positions (Kluth Shuffle)
#                 Returns matrices.
sub randomize_pwm {
    my ($pwm,$verbose)=@_;
    my @pwm;

    if (ref $pwm eq "ARRAY"){
        @pwm = @{$pwm}
    }
    else {
        die "Array parameter is not an array\n";
    }

    
    my @position=(defined $verbose) ? (0..$#pwm) : undef;
    
    
    for ( my $i = $#pwm; $i > 0; --$i ) {
        my $j = int rand( $i + 1 );
        @pwm[ $i, $j ] = @pwm[ $j, $i ];
        if (defined $verbose){
            @position[ $i, $j ] = @position[$j,$i];
        } 
    }
    
    if (defined $verbose){
        print "@position\n";
    #if (main::_query(\@position,\@old_position)){
    #    die "Same matrix was output";
    #}
    }
    
    return @pwm;
}



#Implement scoring function (use threshold and different scoring functions)
#Functions:  Log-likelihood, information content, bits-sub-optimal scoring
#Functions that use scoring are motif masking function, motif identification(R graphing, text output, SVG output )
#Work depreciate function to be used with main functions

###############################################################################
#    Seq_score
#
#   Description: Scoring based upon http://en.wikipedia.org/wiki/Position-specific_scoring_matrix
#
#   Usage: $sequence->Seq_Score($Motif,,"-HTML");
#
#   FLAGS:  -HTML
#   need to fix the $cutoff
#
#   Work:  Depreciate dependency on SEQ class, include SEQ_FREQ dependency
#           Motifset class input:  sequence
#           Returns file with only scores for each motif
#   FIX:  Add option for multiple scoring functions
#   FIX:  Output (verbose, file, variable)
#   FIX:  Threshold (output dependent and function to approximate)
#   FIX:  Automate with R function for graphing motif locations

sub _seq_score {
    #set defaults 
    my %default=("MOTIF"=>undef,
                 "SEQUENCE"=>undef,
                 "SCORING_METHOD"=>"BSO",
                 "THRESHOLD"=>0,
                 "SEQ_FREQ"=>[0.25,0.25,0.25,0.25],
                 );
    my %arg=(%default,@_);
    
    #determine that required variables are supplied
    if (ref $arg{MOTIF} ne "Motif"){
        die "No motif supplied :  $!\n";
    }
    elsif (!defined $arg{SEQUENCE}){
        die "No sequence supplied:  $!\n";
    }
    
    #if passed valid seq_freqency use this for LL and IC else ignore
    if (ref $arg{SEQ_FREQ} eq "HASH"){

    }
    
    my $seq_size=length $arg{SEQUENCE};  #sequence length
    my @motif=@{$arg{MOTIF}->{POSITION_PROBABILITY_MATRIX}};
    my @report;
    
    for (my $i=0;$i<=$seq_size-((scalar @motif)-1);$i++){  #Foreach position along sequence generate a score for both forward Motif and reverse compliment Motif
        my $cumulative_score;
        for (my $j=0;$j<(scalar @motif);$j++){ #Define values and accumulate score foreach position of Motif
            my $position_score;
            my $letter=substr($arg{SEQUENCE},$j+$i,1);
            my $value=  (uc($letter) eq "A") ? ($motif[$j][0]):
                        (uc($letter) eq "C") ? ($motif[$j][1]):
                        (uc($letter) eq "G") ? ($motif[$j][2]):
                        ($motif[$j][3]);
                        
            #start scoring sequence
            #3 different scoring methods  (LL-Log Likelihood, IC-Information Content, BSO- bits-Sub-optimal)
            if ($arg{SCORING_METHOD} eq "LL"){
                $position_score=log($value/$arg{SEQ_FREQ}->{$letter})/log(2);  #compute log likelihood for forward Motif at position
            }
            elsif ($arg{SCORING_METHOD} eq "IC"){
                $position_score=log($value/$arg{SEQ_FREQ}->{$letter})/log(2);  #compute log likelihood for forward Motif at position
                $position_score=$value*$position_score;
            }
            elsif ($arg{SCORING_METHOD} eq "BSO"){
                my $top_value=_max_array(@{$motif[$j]});
                $position_score=(log($value)-log($top_value))/log(2);
            }
            else {
                die "Unsupported scoring method, $arg{SCORING_METHOD}: $!\n";
            }
            $cumulative_score+= $position_score;
        }
        if ($cumulative_score >= $arg{THRESHOLD}){
            #report if above threshold
            push @report, [$i, $i+(scalar @motif),$cumulative_score];
        }
        else {next;}
    }
    return @report;
}

###############################
# mask_motifs
# Description
#
#
#Fix:  Optimize code in for loop better
sub mask_motifs {
    my %default=("MOTIFS"=>undef,
                 "SEQUENCES"=>undef,
                 "SCORING_METHOD" => "BSO",
                 "THRESHOLD"=>0,
                 "SEQ_FREQ"=>[0.25,0.25,0.25,0.25],
                 "STRAND"=>"FORWARD",
                 "OUTPUT"=>undef
                 );
    
    my %arg=(%default,@_);
    
    if (ref $arg{MOTIFS} ne "Motifset"){
        die "mask_motifs function requires Motifset object: $!\n";
    }
    elsif (!defined $arg{SEQUENCES}){
        die "mask_motifs function requires sequence filename: $!\n";
    }
    elsif (!defined $arg{OUTPUT}){
        die "mask_motifs function requires output filename: $!\n";
    }
    #Fix 
    #elsif (($arg{STRAND} ne "FORWARD")||($arg{STRAND} ne "REVERSE")||($arg{STRAND} ne "BOTH")){
    #    die "Strand selection doesn't match:  FORWARD, REVERSE, OR BOTH : $!\n";
    #}
    
    open SEQ, "<$arg{SEQUENCES}";
    open OUTPUT, "> $arg{OUTPUT}";
    my $fasta=new FAlite(\*SEQ);
    while (my $entry=$fasta->nextEntry){
        my $sequence=$entry->seq;
        my $header=$entry->def;
        my $rev_seq=undef;
        my @all_matches;

        #if reverse strand then use only reverse
        if (($arg{STRAND} eq "REVERSE")||($arg{STRAND} eq "BOTH")){
            $rev_seq=reverse $sequence;
            $rev_seq=~tr/ACGTacgt/TGCAtgca/;
        }   


        my $match_function= sub {
            my ($motif,$seq,$scoring_method,$threshold,$seq_freq, $strand)=@_;
            #analyze sequences for motifs
            my @match=_seq_score("MOTIF"=>$motif,
                                   "SEQUENCE"=>$seq,
                                   "SCORING_METHOD"=>$scoring_method,
                                   "THRESHOLD"=>$threshold,
                                   "SEQ_FREQ"=>$seq_freq);
            
            #Corrects the coordinate system for reverse strand
            if ($strand eq "REVERSE"){
                my $seq_size=(length $sequence);
                @match=map ({my @array=@$_;
                               $array[0]=abs($array[0]-$seq_size);
                               $array[1]=abs($array[1]-$seq_size);
                               @array[0,1,2]=@array[1,0,2];
                               $_=\@array;} @match)
            }
            return @match;
        };
        #Find motif matches for given threshold for all motifs supplied
        foreach my $motif (@{$arg{MOTIFS}->{MOTIFS}}){
            my @matches;
            if ($arg{STRAND} eq "FORWARD"){
               @matches=&$match_function($motif,$sequence,$arg{SCORING_METHOD},$arg{THRESHOLD},$arg{SEQ_FREQ},"FORWARD");
               push @all_matches,@matches;
            }
            elsif ($arg{STRAND} eq "REVERSE"){
                @matches=&$match_function($motif,$rev_seq,$arg{SCORING_METHOD},$arg{THRESHOLD},$arg{SEQ_FREQ},"REVERSE");
                push @all_matches,@matches;
            }
            else {
                @matches=&$match_function($motif,$sequence,$arg{SCORING_METHOD},$arg{THRESHOLD},$arg{SEQ_FREQ},"FORWARD");
                push @all_matches,@matches;
                @matches=&$match_function($motif,$rev_seq,$arg{SCORING_METHOD},$arg{THRESHOLD},$arg{SEQ_FREQ},"REVERSE");
                push @all_matches,@matches;
            }
        }
        
        #Consolidate matches (if there are any overlapping borders of motifs), also for next step from 3' end to 5' end deletions
        @all_matches=_organize_masking_regions(\@all_matches);
        
        #delete sequences falling within borders
        my $new_sequence=$sequence;
        foreach my $region (@all_matches){
            if (!defined $region){
                next;
            }
            else{
                my $motif_size=$$region[1] - $$region[0];
                $new_sequence=substr($new_sequence,0,$$region[0]) . "N" x $motif_size . substr($new_sequence,$$region[1],((length $new_sequence)-$$region[1]));
            }
        }
        print OUTPUT "$header". " \(Masked for motifs\)\n";
        print OUTPUT "$new_sequence\n";
        #print "Size of old file: " . length($sequence) . "\n" . "Size of new file: " . length($new_sequence) . "\n";
        #print "$sequence\n$new_sequence"
    }
    close SEQ;
    $|=1;
    close OUTPUT;
}


###################################
# scan_sequence
#
# Description:  Scans sequence file for motifs and reports to standard output
#
#
sub scan_sequence {
    my %default=("MOTIFS"=>undef,
                 "SEQUENCES"=>undef,
                 "SCORING_METHOD" => "BSO",
                 "THRESHOLD"=>0,
                 "SEQ_FREQ"=>[0.25,0.25,0.25,0.25],
                 "STRAND"=>"FORWARD"
                 );
    
    my %arg=(%default,@_);
    
    if (ref $arg{MOTIFS} ne "Motifset"){
        die "mask_motifs function requires Motifset object: $!\n";
    }
    elsif (!defined $arg{SEQUENCES}){
        die "mask_motifs function requires sequence or sequence filename: $!\n";
    }
    #Fix: strand selection
    #elsif (($arg{STRAND}  "FORWARD")||($arg{STRAND} ne "REVERSE")||($arg{STRAND} ne "BOTH")){
    #    die "Strand selection doesn't match:  FORWARD, REVERSE, OR BOTH : $!\n";
    #}
    open SEQ, "<$arg{SEQUENCES}";
    
    my $fasta=new FAlite(\*SEQ);
    while (my $entry=$fasta->nextEntry){
        my $sequence=$entry->seq;
        my $header=$entry->def;
        my $rev_seq=undef;

        #if reverse strand then use only reverse
        if (($arg{STRAND} eq "REVERSE")||($arg{STRAND} eq "BOTH")){
            $rev_seq=reverse $sequence;
            $rev_seq=~tr/ACGTacgt/TGCAtgca/;
        }   


        my $match_function= sub {
            my ($motif,$seq,$scoring_method,$threshold,$seq_freq, $strand)=@_;
            #analyze sequences for motifs
            my @match=_seq_score("MOTIF"=>$motif,
                                   "SEQUENCE"=>$seq,
                                   "SCORING_METHOD"=>$scoring_method,
                                   "THRESHOLD"=>$threshold,
                                   "SEQ_FREQ"=>$seq_freq);
            
            #Corrects the coordinate system for reverse strand
            if ($strand eq "REVERSE"){
                my $seq_size=(length $sequence);
                @match=map ({my @array=@$_;
                               $array[0]=abs($array[0]-$seq_size);
                               $array[1]=abs($array[1]-$seq_size);
                               @array[0,1,2]=@array[1,0,2];
                               push @array, "-";
                               $_=\@array;} @match)
            }
            else {
                @match=map({my @array=@$_;
                            push @array, "+";
                            $_=\@array;} @match)
            }
            
            return @match;
        };
        
        my $print_report=sub{
            my ($arrays, $motif)=@_;
            print "Results for $motif:\n";
            printf "%-10s%-10s%-10s%-10s\n","Start","Finish","Strand","Score";
            print "-" x 40 . "\n";
            if (!defined @$arrays){
                print "NONE\n";
            }
            else{
                foreach my $match (@$arrays){
                    printf "%-10d%-10d%-10s%-8.2f\n",$match->[0]+1,$match->[1],$match->[3],$match->[2];
                }   
            }
        };
        
        
        #Find motif matches for given threshold for all motifs supplied
        foreach my $motif (@{$arg{MOTIFS}->{MOTIFS}}){
            my @matches;
            my @all_matches;
            if ($arg{STRAND} eq "FORWARD"){
               @matches=&$match_function($motif,$sequence,$arg{SCORING_METHOD},$arg{THRESHOLD},$arg{SEQ_FREQ},"FORWARD");
               &$print_report(\@matches,$motif->{ACCESSION});
            }
            elsif ($arg{STRAND} eq "REVERSE"){
                @matches=&$match_function($motif,$rev_seq,$arg{SCORING_METHOD},$arg{THRESHOLD},$arg{SEQ_FREQ},"REVERSE");
                &$print_report(\@matches,$motif->{ACCESSION});
            }
            else {
                @matches=&$match_function($motif,$sequence,$arg{SCORING_METHOD},$arg{THRESHOLD},$arg{SEQ_FREQ},"FORWARD");
                push @all_matches,@matches;
                @matches=&$match_function($motif,$rev_seq,$arg{SCORING_METHOD},$arg{THRESHOLD},$arg{SEQ_FREQ},"REVERSE");
                push @all_matches,@matches;
                @all_matches=sort{${$a}[0] <=> ${$b}[0]} (@all_matches);
                &$print_report(\@all_matches,$motif->{ACCESSION});
            }
        }
    }
    close SEQ;
}







#Takes list of matching regions from each motif and consolidates all overlapping regions into 1 larger region
#Sorts list from 3' right border in decending order
sub _organize_masking_regions {
    my @array=@{shift @_};
    @array=sort {${$b}[1] <=> ${$a}[1]} (@array);
    my @new_array;
    
    #Combine terms and extend borders
    push @new_array, $array[0];
    for (my $i=1;$i<scalar @array;$i++){
        if ($array[$i][1]>$new_array[-1][0]){
            $new_array[-1][0]=$array[$i][0]
        }
        else {
            push @new_array, $array[$i]
        }
    }

    return @new_array;
}



sub scoring_threshold{
    #set defaults
    my %default=("MOTIF"=>undef,
                 "SEQUENCE"=>undef,
                 "SCORING_METHOD"=>"BSO",
                 "SEQ_FREQ"=>[0.25,0.25,0.25,0.25],
                 );
    my %arg=(%default,@_);
    
    #determine that required variables are supplied
    if (ref $arg{MOTIF} ne "Motif"){
        die "No motif supplied :  $!\n";
    }
    elsif (!defined $arg{SEQUENCE}){
        die "No sequence supplied:  $!\n";
    }
    
    #if passed valid seq_freqency use this for LL and IC else ignore
    if (ref $arg{SEQ_FREQ} eq "HASH"){

    }
    
    my $seq_size=length $arg{SEQUENCE};  #sequence length
    my @motif=@{$arg{MOTIF}->{POSITION_PROBABILITY_MATRIX}};
    my @report;
    
    for (my $i=0;$i<=$seq_size-((scalar @motif)-1);$i++){  #Foreach position along sequence generate a score for both forward Motif
        my $cumulative_score;
        for (my $j=0;$j<(scalar @motif);$j++){ #Define values and accumulate score foreach position of Motif
            my $position_score;
            my $letter=substr($arg{SEQUENCE},$j+$i,1);
            my $value=  (uc($letter) eq "A") ? ($motif[$j][0]):
                        (uc($letter) eq "C") ? ($motif[$j][1]):
                        (uc($letter) eq "G") ? ($motif[$j][2]):
                        ($motif[$j][3]);
                        
            #start scoring sequence
            #3 different scoring methods  (LL-Log Likelihood, IC-Information Content, BSO- bits-Sub-optimal)
            if ($arg{SCORING_METHOD} eq "LL"){
                $position_score=log($value/$arg{SEQ_FREQ}->{$letter})/log(2);  #compute log likelihood for forward Motif at position
            }
            elsif ($arg{SCORING_METHOD} eq "IC"){
                $position_score=log($value/$arg{SEQ_FREQ}->{$letter})/log(2);  #compute log likelihood for forward Motif at position
                $position_score=$value*$position_score;
            }
            elsif ($arg{SCORING_METHOD} eq "BSO"){
                my $top_value=_max_array(@{$motif[$j]});
                $position_score=(log($value)-log($top_value))/log(2);
            }
            else {
                die "Unsupported scoring method, $arg{SCORING_METHOD}: $!\n";
            }
            $cumulative_score+= $position_score;
        }
        push @report, [$i, $i+(scalar @motif),$cumulative_score];
    }
    return @report;
}









sub R_seq_score {
    my ($motifset,$sequence,$filename,$seq_freq)=@_;
    my $report=_seq_score($motifset,$sequence,$seq_freq);

    open RFILE, "> $filename.csv";
    my $header="POSITION";
    my @data;
    my $size=length($sequence);

    foreach my $score (@{$report->{SCORES}}){
        my $name=$score->{MOTIF};
        $header.="," . $name ."_LL_FOR";
        $header.="," . $name ."_LL_REV";
        $header.="," . $name ."_IC_FOR";
        $header.="," . $name ."_IC_REV";
        for(my $i=0;$i<$size;$i++){
            if (defined $score->{LL_FOR}->[$i] ){
                $data[$i].="," . $score->{LL_FOR}->[$i];
                $data[$i].="," . $score->{LL_REV}->[$i];
                $data[$i].="," . $score->{IC_FOR}->[$i];
                $data[$i].="," . $score->{IC_REV}->[$i];
            }
            else{
                $data[$i].="," ;
                $data[$i].="," ;
                $data[$i].="," ;
                $data[$i].="," ;
            }

        }
    }

    for (my $i=1;$i<$size;$i++){
        $data[$i]=$i . $data[$i];
    }

    print RFILE "$header\n";
    my $data=join("\n",@data);
    print RFILE $data;

    close RFILE;
    $|=1
}



############################################
#
#   _query()
#
#   Compares the variable to the search term
#   $searchterm can be scalar, hash, array, reference
#
#
#
sub _query{
    my ($variable, $searchterm)=@_;
    my $return_value=0;

    if (!ref $searchterm){
        if (ref $variable eq "HASH"){
            foreach my $key (keys %$variable){
                if ($key eq $searchterm){
                    $return_value=1;
                    last;
                }
                elsif (_query($variable->{$key},$searchterm)){
                    $return_value=1;
                    last;
                }
                else {
                    next;
                }
            }

        }
        elsif (ref $variable eq "ARRAY"){
            foreach my $array (@$variable){
                if (_query($array,$searchterm)){
                    $return_value=1;
                    last;
                }
            }

        }
        elsif (ref $variable eq "SCALAR"){
            if ($$variable eq $searchterm){
                $return_value=1;
            }
        }
        elsif (ref $variable eq "REF"){
            if(_query($$variable,$searchterm)){
                $return_value=1;
            }
        }
        elsif ($variable=~/$searchterm/g){
           $return_value=1;
        }
    }

    elsif (ref $searchterm eq "ARRAY"){
        if (ref $variable eq "HASH"){
            foreach my $key (keys %$variable){
                if ($key eq $searchterm){
                    $return_value=1;
                    last;
                }
                elsif (_query($variable->{$key},$searchterm)){
                    $return_value=1;
                    last;
                }
                else {
                    next;
                }
            }

        }
        elsif (ref $variable eq "ARRAY"){
            my @array_search=@$searchterm;
            my @array_var=@$variable;

            if (scalar @array_search == scalar @array_var){
                my $str1=join('',@array_search);
                my $str2=join('',@array_var);
                if ($str1 eq $str2){
                    $return_value=1;
                }
            }
        }
        elsif (ref $variable eq "REF"){
            if (_query($$variable,$searchterm)){
                $return_value=1;
            }
        }
    }

    elsif (ref $searchterm eq "HASH"){
        if (ref $variable eq "HASH"){
            my %searchterm=%$searchterm;
            my %variable=%$variable;
            my $search_keys=join('',sort keys %searchterm);
            my $variable_keys=join('',sort keys %variable);
            my $search_values=join('',sort values %searchterm);
            my $variable_values=join('',sort values %variable);
            if ((scalar $search_keys == scalar $variable_keys) && (scalar $search_values == scalar $variable_values)){
                if (($variable_keys eq $search_keys)&&($variable_values eq $search_values)){
                    $return_value=1;
                }
            }
        }
        elsif (ref $variable eq "ARRAY"){
            foreach (@$variable){
                if (_query($_,$searchterm)){
                    $return_value=1;
                    last;
                }
            }
        }
        elsif (ref $variable eq "REF"){
            if (_query($$variable,$searchterm)){
                $return_value=1;
            }
        }
    }

    elsif (ref $searchterm eq "REF"){
        if (_query($variable,$$searchterm)){
            $return_value=1;
        }
    }

  return $return_value;
}

sub svg_template {
    my @svg=('<?xml version="1.0" standalone="no"?>',
             '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">',
             '<svg width="$length" height="150" version="1.1" xmlns="http://www.w3.org/2000/svg">',
             '<text transform="matrix(1 0 0 $transformation $horizontal_position $vertical_position)" fill="$letter_color" font-family="Helvetica" font-size="120" text-anchor="middle">$letter</text>',
             '<!-- Setup graph and add axis labels -->',
             '<path d="M 30,30 V 120 H $graphlength" fill="none" stroke="black"/>',
             '<path d="M 25,75 H 30" stroke="black"/>',
             '<path d="M 25,30 H 30" stroke="black"/>',
             '<path d="M 25,120 H 30" stroke="black"/>',
             '<!-- Y axis labels and tick marks -->',
             '<text transform="matrix(1 0 0 1 18 125)" font-family="Helvetica-Bold" font-size="12">0</text>',
             '<text transform="matrix(1 0 0 1 12 80)" font-family="Helvetica-Bold" font-size="12">50</text>',
             '<text transform="matrix(1 0 0 1 6 35)" font-family="Helvetica-Bold" font-size="12">100</text>',
             '<text transform="matrix(0 -1 1 0 10 110)" font-family="Helvetica-Bold" font-size="12">Percentage</text>',
             '<!-- X axis labels -->',
             '<text transform="matrix(1 0 0 1 25 135)" font-family="Helvetica-Bold" font-size="12">5</text>',
             '<text transform="matrix(1 0 0 1 $graph_length 135)" font-family="Helvetica-Bold" font-size="12">3</text>',
             '<!-- X axis labels-->',
             '<text transform="matrix(1 0 0 1 $position_segment 135)" font-family="Helvetica-Bold" font-size="12" text-anchor="middle">$position</text>',
             '</svg>');
    return \@svg;
}

sub nm_template {
    my @nm=('<motifset>',
            '   <prop>',
            '       <key>$key</key>',
            '       <value>$value</value>',
            '   </prop>',
            '   <motif>',
            '       <name>$motif</name>',
            '       <weightmatrix alphabet="DNA" columns="$columns">',
            '           <column pos="$column">',
            '               <weight symbol="adenine">$a_value</weight>',
            '               <weight symbol="cytosine">$c_value</weight>',
            '               <weight symbol="guanine">$g_value</weight>',
            '               <weight symbol="thymine">$t_value</weight>',
            '           </column>',
            '       </weightmatrix>',
            '       <threshold>$threshold</threshold>',
            '   </motif>',
            '</motifset>');
    return \@nm;
}


sub html_index {
    my @index_html=('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Frameset//EN" "http://www.w3.org/TR/html4/frameset.dtd">',
              '<html>',
              '<head>',
              '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1">',
              '<title>Motif Viewer</title>',
              '</head>',
              '<frameset rows="150,150,150,150,150" cols="200,*" framespacing="0" frameborder="yes" border="2">',
              '<frame src="Menu1.html">',
              '<frame src="Summary.html" name="Motif1" scrolling="No" noresize="noresize" id="Motif1" title="Motif1" >',
              '<frame src="Menu2.html">',
              '<frame src="Summary.html" name="Motif2" scrolling="No" noresize="noresize" id="Motif2" title="Motif2" >',
              '<frame src="Menu3.html">',
              '<frame src="Summary.html" name="Motif3" scrolling="No" noresize="noresize" id="Motif3" title="Motif3" >',
              '<frame src="Menu4.html">',
              '<frame src="Summary.html" name="Motif4" scrolling="No" noresize="noresize" id="Motif4" title="Motif4" >',
              '<frame src="Menu5.html">',
              '<frame src="Summary.html" name="Motif5" scrolling="No" noresize="noresize" id="Motif5" title="Motif5" >',
              '</frameset>',
              '<noframes><body>',
              '</body>',
              '</noframes></html>'
              );
    return \@index_html;  
}

sub html_menu{
    my @menu_html=('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">',
                   '<html xmlns="http://www.w3.org/1999/xhtml">',
                   '<head>',
                   '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1" />',
                   '<title>Index Menu</title>',
                   '<SCRIPT TYPE="text/javascript">',
                   '<!--',
                   'function dropdown(mySel)',
                   '{',
                   'var myWin, myVal;',
                   'myVal = mySel.options[mySel.selectedIndex].value;',
                   'if(myVal)',
                    '    {',
                    '    if(mySel.form.target)myWin = parent[mySel.form.target];',
                    '    else myWin = window;',
                    '    if (! myWin) return true;',
                    '    myWin.location = myVal;',
                    '    }',
                    'return false;',
                    '}',
                    '//-->',
                    '</SCRIPT>',
                    '</head>',
                    '<body>',
                    '<FORM',
                    '    ACTION="Summary.html"',
                    '    METHOD=POST onSubmit="return dropdown(this.gourl)"',
                    '    TARGET={$Motif}',
                    '    >',
                    '<SELECT NAME="gourl">',
                    '</SELECT>',
                    '<INPUT TYPE=SUBMIT VALUE="Go">',
                    '</FORM>',
                    '</body>',
                    '</html>',
                   );
    return \@menu_html;
}

sub html_summary{
    my @summary_html=('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 TRANSITIONAL//EN">',
                      '<html>',
                      '     <head>',
                      '     <title>{$JOB}</title>',
                      '     </head>',
                      '     <body>');
    return \@summary_html;
}

#Returns the maximum value of an array
sub _max_array{
    my @array=@_;
    my $max=$array[0];
    foreach (@array){
        if ($max<$_){
            $max=$_;
        }
    }
    return $max;
}

1;

__END__


=head1 NAME
Motif  - Motif module

=head1 SYNOPSIS

    use Motif;

    my $Demo_Motifset=Motifset->_new_Motifset();
    my $meme=Motifset->_new_Motifset();
    $Demo_Motifset->Transfac_import("Matrix.dat");
    $Demo_Motifset->NestedMica_import("Nested\ Mica\ Motifs.xms");
    $meme->Meme_import("Meme\ Output\ DNA.txt");
    $Demo_Motifset->Weeder_import("G1S.cycle.fasta.wee");
    $Demo_Motifset->PrintSummary(-a);
    $Demo_Motifset->combine_Motifset($meme);
    $Demo_Motifset->save_data("MyWeeder.txt");
    my $retrieved_data=Motifset->retrieve_data("MyWeeder.txt");
    $Demo_Motifset->Splitstree("NM","Job2");
    $Demo_Motifset->Weblogo_Summary("Job2");
    my $sequences=SEQUENCE_SET->new();
    $sequences->Fasta_Import("ADFP.txt");
    my $new_Motifset=Motifset->_new_Motifset();
    $new_Motifset->subset($Demo_Motifset);
    $new_Motifset->delete_Motif();

    my $seq=SEQUENCE->new();
    $seq->simple_shuffle($sequences,5,5);
    $seq->randomize($sequences);
    my $scoreset=$Demo_Motifset->Seq_Score($sequences);

    $sequence->generate_sequence($sequence,1700);
    $sequence->embed_Motif($sequence,$Demo_Motifset,100,3);
    $sequence->embed_Motif($sequence,$Demo_Motifset,NA,3);
    $sequence->embed_Motif($sequence,$Demo_Motifset,100,NA);
    $sequence->embed_Motif($sequence,$Demo_Motifset);
    my $seq=SEQUENCE_SET->new();
    $seq->subset($sequences);
    $seq->delete_sequence();


=head1 Description

Motif is a package for analysis and storing Motif data from (Transfac, Nested_Mica, Weeder, Meme).
Importing of data is stored in a common internal format as objects; described in the next section.
There are many different Motif finders each with there own strengths.   However, there wasn't a
utilitie to organize, analyse, and combine data from different Motif finders under a common format.

Motif contains many of tools for storing, analyzing Motifs into one compact package.

=head2 Objects

There are 4 different objects
Motif :  Each individual Motif contains
= over 4
ACCESSION = unique identifier for each Motif.  Transfac accessions are used as accession.  However, all other Motifs use accessions that are generated in the format "NMplott82807N00"  Meaning the Motif was imported from NestedMica by user "plott" on Aug 28,2007.
IDENTIFIER = Motif name from file
Motif_JOB = stores Meme Job number that Motif comes from
DATE_AUTHOR_UPDATE_DATES
COPYRIGHTNAME_OF_BINDING_SITE
SHORT_FACTOR_DESCRIPTION
LIST_OF_LINKED_FACTOR_ENTRIES
STATISTICAL_BASIS
BINDING_SITES_UNDERLYING_MATRIX
COMMENTS
MEDLINE_ID
REFERENCE_#
REFERENCE_AUTHORS
REFERENCE_TITLE
REFERENCE_DATA
Motif_SEQUENCES
SEQUENCES
Motif_METHOD
COMMAND_LINE
DATAFILE
ALPHABET
REGULAR_EXPRESSION
CONSENSUS_SEQUENCE
VERSION
THRESHOLD
ORGANISM
POSITION_SCORING_MATRIX
POSITION_PROBABILITY_MATRIX
BACKGROUND_FREQUENCY

=back

Motifset= set of Motifs, can also contain Distance_Score (which is a measurement between Motifs)

SEQUENCE= contains to following:
=over 4
DATABASE: database the sequence was from - information gotten from FASTA header;  if not known it is set to "Unknown"
ACCESSION: Accession from FASTA
LOCUS:  Locus information from header
FORWARD: Forward sequence
REVERSE_COMP: Reverse compliment sequence
SEQ_FREQ:  Sequence mono,di,tri nucleotide frequencies

=back

SEQUENCE_SET:  Set of sequence objects.  Contains:
=over 4
SEQUENCE:
SEQ_FREQ:Sequence mono,di,tri nucleotide frequencies for the complete set
=back


SEQUENCE_FREQUENCY: contains frequency table; by default mono,di,tri nucleotides are measured


Distance_Score: distance scores between Motifs.  Contains:
=over 4
METHOD = method used to create distance score
DISTANCE_MATRIX = matrix of distances
LABEL = Labels for table, contains list of Motifs used to create distance matrix
=back


=head1 Package Summary



=head2 Package: Motif Function Summary
=over 4
=item * new_Motif() : Returns new empty Motif object;
=item *
=back

=head2 Package: Motifset Function Summary

=over 4
=item * new_Motifset() :  Returns empty Motifset object
=item * Transfac_import("Filename") or Transfac_import()  :  Import transfac matrix data file.   If filename not defined, query user for filename. Returns Motifset of Motifs imported from filename
=item * NestedMica_import("Filename") or NestedMica_import(): Import NestedMica xms data file.  Returns Motifset of Motifs imported from  file
=item * Meme_import("Filename") or Meme_import():   Import Meme text data file.  Returns Motifset of Motifs imported from file (Doesn't handle HTML format)
=item * Weeder_import("Filename") or Weeder_import():  Import Weeder data file.  Returns Motifset of Motifs imported from file.
=item * PrintSummary(-v) or PrintSummary(-a):  Flag -v prints verbose summary (everything stored in Motif); -s prints short summary( Accession, Idenifier, PPM and PSSM)
=item * combine_Motifset($Motifset):  Combines two Motifsets and returns the combined set
=item * save_data($filename) : Stores Motifset to file
=item * retrieve_data($filename) : Retrieves Motifset from file
=item * import_Motif() : not fully implemented
=item * _generate_consensus: not fully implemented
=item * Score_Motifset("Method"):  Methods: NM=nested mica, KL=Kullback-Leibler Distance  other scoring methods not fully implemented
=item * Splitstree("Method","DirectoryName"): creates Nexus output file for Splitstree
=item * Weblogo_Summary("DirectoryName")  : creates sequence logo with frequency rather than IC, Creates index.html for easy viewing and comparing of logos
=item * subset(Motifset) : takes a Motifset and returns a subset of the that Motifset
=item * delete_Motif() :  Deletes Motif from the Motifset
=back

=head2 Package: PSM_MATRIX
=over 4
=item * CONVERT(): Converts PSM to PPM.  Returns PPM
=back

=head2 Package: PPM_MATRIX
=over 4
=item * new_matrix():  creates new Position Probability Matrix object
=item * NM_score();  used by Score_Motifset to measure distances between Motifs
=item * NM_functon();  generates score for NM_score
=back

=head2 Package:  SEQ_FREQ
=over 4
=item * new() : creates new SEQ_FREQ object
=item * frequency(SEQUENCE) or frequency(SEQUENCE_SET): returns frequency table for sequence or sequence set
=back

=head2 Package: Distance_Score
=over 4
=item * new(DISTANCE_MATRIX,LABEL,METHOD) : returns Distance_Score object;  called by Motifset->Score_Motifset
=back

=head2 Package: SEQUENCE_SET;
=over 4
=item * new() : returns empty SEQUENCE_SET object
=item * Fasta_Import("filename"): imports Fasta sequences into SEQUENCE_SET
=item * save_fasta("filename"): saves the SEQUENCE_SET to fasta file
=item * save_sequence_set("filename"):  saves SEQUENCE_SET in internal format
=item * retrieve_sequence_set("filename"):  retrieves SEQUENCE_SET from file
=item * subset(SEQUENCE_SET) : returns a subset of the sequence set
=item * delete_sequence() : deletes user selected sequence from sequence set
=back

=head2 Package: SEQUENCE
=over 4
=item * new(): creates new empty SEQUENCE object
=item * reverseComplement(): creates reverse compliment of sequence in sequence object
=item * simple_shuffle($sequence,$n_shuffles,$word_size): takes sequence object, number of shuffles, and word size to use in shuffle and returns seqence object that is shuffled like deck of cards
=item * randomize($sequence): takes a sequence object and returns a sequence object that is randomized using Knuth shuffle
=item * save_fasta("filename"): saves Sequence object as Fasta file
=item * save_sequence("filename"): saves Sequence object to file
=item * retrieve_sequence("filename"): retrieves sequence object from file
=item * generate_sequence($sequence,size): generates sequence using single nucleotide frequency table from sequence or sequence set  #need to impliment background frequency from Motifset/Motif.....
=item * extract_sequence($sequence_set): takes a sequence_set object and returns a sequence after querying the user as to the accession ot extract
=item * Seq_Score(Motif/Motifset,cutoff,flag): if given sequence will use background info from sequence to generate sequence scores.   Returns information and HTML summary of matching Motif positions.
=item * _HTML_report(   ):  used by Seq_score to report scoring in Summary.html
=item * embed_Motif($sequence,$Motif,$spacing,$n): takes a sequence object, a Motif( if given Motifset will query user), spacing size, number of Motifs; and returns sequence object with Motif embedded.  Uses random number generator to generate each embedded Motif from PPM then imports that Motif at one position
=back

=head2 Package: main
=over 4
=item * _query_user();  If handed a set will query user for selection and return selection
=back

=head1 Function Descriptions and Examples

=head2 Motif Package

=head3 new_Motif

=over4
INPUT:  NONE
OUTPUT: Empty Motif Object

EXAMPLE:  my $Motif->Motif->new_Motif();

=back

=head3 clone

=over4
INPUT:  Motif
OUTPUT: Empty Motif Object

EXAMPLE:  my $cloned_Motif->clone($Motif);
$cloned_Motif is a clone of $Motif
=back


=head2 Motifset

=head3 _new_Motifset()

=over4
INPUT:  NONE
OUTPUT: Empty Motifset Object

EXAMPLE:  my $Motifset->Motifset->_new_Motifset();
=back

=head3  PrintSummary

=over4
INPUT: Motifset,flag
FLAGS: -S,-V, none
OUTPUT: Output to STDOUT summary of Motifset information

EXAMPLE:

my $Motifset->PrintSummary();
=over5
-----------------------------------------------
Motif: MEplott9177N00
-----------------------------------------------
ACCESSION: MEplott9177N00
Import File: Mfn.txt

-------------------------
Position Probability Matrix:
-------------------------
A    C    G    T
0.00 0.25 0.00 0.75
0.12 0.50 0.38 0.00
0.25 0.00 0.25 0.50
0.12 0.25 0.12 0.50
0.00 0.38 0.62 0.00
0.00 0.38 0.00 0.62
0.25 0.25 0.38 0.12
0.25 0.25 0.12 0.38
0.25 0.12 0.38 0.25
0.25 0.00 0.00 0.75
0.00 0.00 1.00 0.00
0.00 0.00 0.50 0.50
0.00 0.25 0.62 0.12
0.62 0.12 0.00 0.25
0.00 0.00 0.25 0.75
0.00 0.00 0.50 0.50
0.00 0.00 0.75 0.25
0.00 0.00 0.75 0.25
0.00 0.50 0.50 0.00
0.00 0.75 0.25 0.00
0.25 0.62 0.00 0.12
0.38 0.38 0.12 0.12
0.75 0.00 0.25 0.00
0.50 0.38 0.00 0.12
0.38 0.50 0.12 0.00
0.25 0.00 0.38 0.38
0.00 1.00 0.00 0.00
0.25 0.38 0.00 0.38
0.50 0.25 0.25 0.00
1.00 0.00 0.00 0.00
0.25 0.00 0.38 0.38
0.00 0.12 0.00 0.88
0.00 0.38 0.62 0.00
0.25 0.62 0.00 0.12
0.25 0.25 0.12 0.38
0.00 0.75 0.00 0.25
0.38 0.50 0.00 0.12
0.00 1.00 0.00 0.00
0.00 0.00 0.00 1.00
0.25 0.38 0.25 0.12
0.50 0.38 0.00 0.12
0.00 0.00 0.25 0.75
0.00 0.12 0.62 0.25
0.38 0.12 0.25 0.25
0.38 0.00 0.25 0.38
0.25 0.00 0.38 0.38
0.25 0.00 0.75 0.00
0.62 0.38 0.00 0.00
0.12 0.00 0.50 0.38
0.00 0.12 0.88 0.00
.......
=back

my $Motifset->PrintSummary("-S");
=over4
-----------------------------------------------
Motif: MEplott9177N00
-----------------------------------------------
Motif: MEplott9177N01
-----------------------------------------------
Motif: MEplott9177N02
=back

my $Motifset->PrintSummary("-V");
=over 4
-----------------------------------------------
Motif: MEplott9177N00
-----------------------------------------------
ACCESSION: MEplott9177N00
Import File: Mfn.txt

-------------------------
Position Probability Matrix:
-------------------------
A    C    G    T
0.00 0.25 0.00 0.75
0.12 0.50 0.38 0.00
0.25 0.00 0.25 0.50
0.12 0.25 0.12 0.50
0.00 0.38 0.62 0.00
0.00 0.38 0.00 0.62
0.25 0.25 0.38 0.12
0.25 0.25 0.12 0.38
0.25 0.12 0.38 0.25
0.25 0.00 0.00 0.75
0.00 0.00 1.00 0.00
0.00 0.00 0.50 0.50
0.00 0.25 0.62 0.12
0.62 0.12 0.00 0.25
0.00 0.00 0.25 0.75
0.00 0.00 0.50 0.50
0.00 0.00 0.75 0.25
0.00 0.00 0.75 0.25
0.00 0.50 0.50 0.00
0.00 0.75 0.25 0.00
0.25 0.62 0.00 0.12
0.38 0.38 0.12 0.12
0.75 0.00 0.25 0.00
0.50 0.38 0.00 0.12
0.38 0.50 0.12 0.00
0.25 0.00 0.38 0.38
0.00 1.00 0.00 0.00
0.25 0.38 0.00 0.38
0.50 0.25 0.25 0.00
1.00 0.00 0.00 0.00
0.25 0.00 0.38 0.38
0.00 0.12 0.00 0.88
0.00 0.38 0.62 0.00
0.25 0.62 0.00 0.12
0.25 0.25 0.12 0.38
0.00 0.75 0.00 0.25
0.38 0.50 0.00 0.12
0.00 1.00 0.00 0.00
0.00 0.00 0.00 1.00
0.25 0.38 0.25 0.12
0.50 0.38 0.00 0.12
0.00 0.00 0.25 0.75
0.00 0.12 0.62 0.25
0.38 0.12 0.25 0.25
0.38 0.00 0.25 0.38
0.25 0.00 0.38 0.38
0.25 0.00 0.75 0.00
0.62 0.38 0.00 0.00
0.12 0.00 0.50 0.38
0.00 0.12 0.88 0.00

-------------------------
Motif_JOB
-------------------------
20181
-------------------------
SEQUENCES
-------------------------
gi|41350328|ref|NM_01487
gi|45269136|ref|NM_03354
-------------------------
COMMAND_LINE
-------------------------
meme /home/meme/meme354/LOGS/meme.20181.data -dna -mod tcm -nMotifs 3 -minw 6 -maxw 50 -evt 1e100 -time 7200 -maxsize 60000 -nostatus -maxiter 20 -dir /home/meme/meme354
-------------------------
DATAFILE
-------------------------
 pasted_sequences
-------------------------
ALPHABET
-------------------------
 ACGT
-------------------------
REGULAR_EXPRESSION
-------------------------
[TC][CG][TAG][TC][GC][TC][GAC][TAC][GAT][TA]G[GT][GC][AT][TG][GT][GT][GT][CG][CG][CA][AC][AG][AC][CA][GTA]C[CTA][ACG]A[GTA]T[GC][CA][TAC][CT][CA]CT[CAG][AC][TG][GT][AGT][ATG][GTA][GA][AC][GT]G

-------------------------
VERSION
-------------------------
version 3.5.4

=back
=back

=head3  delete_Motif

=over4
INPUT: Motifset, (OPTIONAL: ACCESSION NUMBER)
OUTPUT: Motifset with the user defined Motif deleted from Motifset

EXAMPLE:

$new_Motifset->delete_Motif();  #will query user for the Motif to delete

or

$new_Motifset->delete_Motif(MEplott9177N02);  #will delete MEplott9177N02 Motif from the Motifset


=back

=head3  combine_Motifset()

=over4
INPUT: Motifset, Motifset
OUTPUT: Motifset

EXAMPLE:

$Motifset->combine_Motifset($second_Motifset);

Combines $second_Motifset with $Motifset, so that $Motifset now has all the combined Motifs
=back

=head3  import_Motif

=over4
INPUT: Motifset, Motif
OUTPUT: Motifset with Motif imported

EXAMPLE:

$Motifset->import_Motif($Motif)
=back

=head3  subset_Motifset

=over4
INPUT: Motifset, (OPTIONAL: array of numbers corresponding to index)
OUTPUT: Motifset containing the selected subset

EXAMPLE:

$subset_Motifset->subset_Motifset($Motifset);

1) MEplott9177N00
2) MEplott9177N01
3) MEplott9177N02
Enter the Motifs that you would like to create subset( Example: 1,2,3,4):1,2

Now, $subset contains MEplott9177N00 and MEplott9177N01

If the index values are known,
@Motif_index=(1,2);
$subset_Motifset->subset_Motifset($Motifset,\@Motif_index);

Now, $subset contains MEplott9177N00 and MEplott9177N01

=back


=head3  save_data

=over4
INPUT: "$filename"
OUTPUT: Saves data to filename. Confirmation to STDOUT.

EXAMPLE:

$Motifset->save_data("Motifset1.PCL");
STDOUT: Motifset saved to Motifset1.PCL
=back

=head3  retrieve_data

=over4
INPUT: "$filename"
OUTPUT: Retrieves data from filename. Confirmation to STDOUT.

EXAMPLE:

$Motifset->retrieve_data("Motifset1.PCL");
STDOUT: Motifset retrieved from Motifset1.PCL
=back


=head3  Score_Motifset
Scores all the Motifs in the Motifset against each other to deteremine the distance between Motifs.   The distances are stored in the Motifset as a Distance_Score matrix

=over4
INPUT: Method for scoring
OUTPUT: Creates distance score matrix in Motifset

Scoring Methods:
PCC:  Pearson correllation coefficient
ALLR: Average Log-Likelihood ratio
SN: Sandelin/Wasserman
KL: Kullback-Leibler Distance
NM: Nested Mica distance
NMKL: Kullback-Leibler/Nested Mica Hybrid


EXAMPLE:








































1;