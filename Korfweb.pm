###############################################################################
# Korfweb - standardized web applications in the Korf lab
###############################################################################
package Korfweb;
use strict;
use warnings;
use CGI qw(:standard);
use File::Temp qw(tempfile);

our $STYLE = "\
body {
        color: black; background: white;
        font-family: Verdana,Arial,Trebuchet MS, sans-serif;
}
h1,h2,h3,h4,h5,h6 { 
        font-family: Trebuchet MS, Verdana;
}

pre { font-family: monospace; }
 
p { text-align: justify; }
";

my $CONTACT = 'For bug reports, feature requests, or general feedback, please email <a href="mailto:korflab@ucdavis.edu">korflab@ucdavis.edu</a>. For more information, see our website at <a href=http://korflab.ucdavis.edu>korflab.ucdavis.edu</a>.';

sub new {
	my ($class, %p) = @_;
	my $self = bless {};
	
	die unless defined $p{title};
	die unless defined $p{element};
	$self->{cgi}     = new CGI;
	$self->{title}   = $p{title};
	$self->{info}    = defined $p{info} ? $p{info} : "";
	$self->{refs}    = defined $p{refs} ? $p{refs} : "";
	$self->{element} = $p{element};
	
	return $self;
}

sub display {
	my ($self) = @_;
	
	my $c = $self->{cgi};
	print $c->header;
	print $c->start_html(-title => $self->{title},-style => {-code => $STYLE});

	print "<CENTER>\n";
        print "<a href=\"/\"><IMG alt=\"Korflab\" border=0 src=\"/images/korflab_logo1_large.jpeg\"></A>\n";
        print "<br>\n";
        print "<TABLE bgColor=#0000cd border=0 cellPadding=1 cellSpacing=0>\n";
        print "<TR align=middle>\n";
        print "<TD width=110><A href=\"/about.html\"><FONT color=#ffffff>About</FONT></A></TD>\n";
    	print "<TD width=110><A href=\"/people.html\"><FONT color=#ffffff>People</FONT></A></TD>\n";
    	print "<TD width=110><A href=\"/research.html\"><FONT color=#ffffff>Research</FONT></A></TD>\n";
    	print "<TD width=110><A href=\"/publications.html\"><FONT color=#ffffff>Publications</FONT></A></TD>\n";
    	print "<TD width=110><A href=\"/software.html\"><FONT color=#ffffff>Software & Data</FONT></A></TD>\n";
    	print "<TD width=110><A href=\"/Unix_and_Perl/index.html\"><FONT color=#ffffff>Unix & Perl course</FONT></A></TD>\n";
        print "</TR>\n";
	  	print "<TR bgcolor=#191970 align middle>\n";
	    print "<TD colspan=7 align=center><font color=#ffffff size=\"-1\"><!--#config timefmt=\"%A %B %d, %Y\" --><!--#echo var=\"DATE_LOCAL\" --></font></TD>\n";
		print "</TR>\n";
        print "</TABLE>\n";
        print "</CENTER>\n";

	print $c->h1($self->{title});
	print $self->{info}, $c->p;
	
	# reference section
	if ($self->{refs}) {
		print $c->h3('References');
		print "<ul>";
		foreach my $string (@{$self->{refs}}) {print "<li>", $string, "<br>\n"}
		print "</ul>";
	}
	
	# form section
	print $c->h3('Input');
	print $c->start_multipart_form;
	foreach my $element (@{$self->{element}}) {
		my $type = shift @$element;
		if ($type eq 'fasta') {
			my ($name, $string) = @$element;
			print $string, $c->br;
			print $c->textarea(
				-name => $name,
				-rows => 5,
				-columns => 80), $c->p;
		} elsif ($type eq 'scrollbox') {
			my $name = shift @$element;
			my $text = shift @$element;
			my @value;
			my %label;
			foreach my $item (@$element) {
				my ($key, $string) = (@$item);
				push @value, $key;
				$label{$key} = $string;
			}
			my $size = @value > 10 ? 10 : @value;
			print $text, $c->br,
				$c->scrolling_list(
					-name => $name,
					'-values' => \@value,
					-labels => \%label,
					-size => $size,
					-multiple => 'false'),
				$c->p;
		} elsif ($type eq 'textbox') {
			my ($name, $text, $size, $default) = @$element;
			print $text, " ",
				$c->textfield(
					-name => $name,
					-size => $size,
					-default => $default),
				$c->br;
		} elsif ($type eq 'text') {
			print @$element;
			print $c->br;
		} elsif ($type eq 'radio_group'){
			my $name = shift @$element;
			my @value;
			foreach my $item (@$element) {
				push @value, $item;
			}
			print $text, $c->br,
				$c->radio_group(
					-name => $name,
					'-values' => \@value,),
				$c->p;
		}
	}
	print $c->p;
	print $c->submit, $c->defaults('Clear Form');

	# feedback section
	print $c->hr;
	print $c->h3('Feedback');
	print $CONTACT;
		
	print $c->end_form;
	return $c->Vars;

}

sub is_non_negative_int {
	my ($val) = @_;
	if    (not defined $val) {return 0}
	elsif ($val < 0)         {return 0}
	elsif ($val != int $val) {return 0}
	else                     {return 1}
}

sub is_fasta {
	my ($text) = @_;
	my @line = split(/\n/, $text);
	if ($line[0] !~ />\S+/) {
		return 0;
	}
	for (my $i = 0; $i < @line; $i++) {
		if ($line[$i] =~ /^>/) {
			return 0 if $line[$i+1] !~ /^\w+$/;
		}
	}
	return 1;
}

1;

__END__

