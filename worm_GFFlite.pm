# this perl module is for creating subset files from worm chromosome GFF files
# grabs records of genes, cds, intron, exon, protein_coding_transcripts
# usage: &GFFlite::gff_reader_alt(directory, gff_file);
#
# Written by Yi Zhang, July/August 2009

package worm_GFFlite; 
use strict; 
use warnings;

sub gff_reader_alt {
	
	my ($dir, $gff_file) = @_;
	my $dir_file = "$dir" . "$gff_file";
	my $output_filename = "$dir"."subfile_" . "$gff_file";
	my %CDS_to_status = ();
	my @lines; #stores all the records that we're interested in
	
	open (IN, "<$dir_file") or die "GFFlite: Can't open $dir_file\n";
	
	while (my $line = <IN>) {
		next if ($line =~ /^#/);
		chomp $line;
		my ($seqname, $source, $feature, $start, $stop, $score, $strand, $frame, $comments) = split (/\t/, $line);
		if ($source eq 'curated' and $feature eq 'CDS') {
			$comments =~ m/CDS \"(.*?)\".*/; my $cdsid = $1;
			$comments =~ /Gene \"(WBGene\d{8})\".*/; my $gid = $1;	
			$comments =~ /Status \"(\S+)\".*/; my $status = $1;
			$CDS_to_status{$cdsid} = $status;
			my $revised_comments = "$gid ; "."CDS_STATUS: $status";
			push (@lines, "$gff_file\tCDS\t$cdsid\t$start\t$stop\t$strand\t$revised_comments");
		}
		if ($source eq 'curated' and $feature eq 'intron') {
			my ($cdsid) = $comments =~ m/CDS \"(.*?)\".*/; 
			push (@lines, "$gff_file\tintron\t$cdsid\t$start\t$stop\t$strand");
		}
		if ($source eq 'curated' and $feature eq 'exon') {
			my ($cdsid) = $comments =~ m/CDS \"(.*?)\".*/; 
			push (@lines, "$gff_file\texon\t$cdsid\t$start\t$stop\t$strand");		
		}
		if ($source eq 'gene' and $feature eq 'gene') {
			my ($gid) = $comments =~ /Gene \"(WBGene\d{8})\".*/; 	
			push (@lines, "$gff_file\tgene\t$gid\t$start\t$stop\t$strand");
		}
		if ($source eq 'Coding_transcript' and $feature eq 'protein_coding_primary_transcript') {
			my ($transcript) = $comments =~ /Transcript \"(.*?)\".*/;
			$transcript =~ /(\S+\.\S+).*?/; my $cdsid = $1;
			push (@lines, "$gff_file\ttranscript\t$transcript\t$start\t$stop\t$strand\t$cdsid");
		}
	}
	close (IN); 
	
	# adds the status to the last column of exon, intron, transcript records
	open (OUT, ">$output_filename") or die "GFFlite: Can't open $output_filename to write\n";
	foreach my $line (@lines) {
		my ($chromosome, $feature, $featureid, $start, $stop, $strand, $comment) = split(/\t/, $line);
		if ($feature eq 'intron' or $feature eq 'exon') {
			$line .= "\tCDS_STATUS: $CDS_to_status{$featureid}";
		}
		elsif ($feature eq 'transcript') {
			$line .= " ; CDS_STATUS: $CDS_to_status{$comment}";
		}
		print OUT "$line\n";
	}
	close (OUT);
	return $output_filename; #returns the name and exact location of the created file
}

1;

__END__

