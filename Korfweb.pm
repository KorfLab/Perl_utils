###############################################################################
# Korfweb - standardized web applications in the Korf lab
###############################################################################
package Korfweb;
use strict;
use warnings;
use CGI qw(:standard);
use File::Temp qw(tempfile);

our $STYLE = "\
H1 {
	color: red;
}
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
	print $c->h1($self->{title});
	print $self->{info}, $c->p;
	
	# reference section
	if ($self->{refs}) {
		print $c->h3('References');
		print "<ul>";
		foreach my $string (@{$self->{refs}}) {print "<li>", $string, "<br>\n"}
		print "</ul>";
	}
	
	# feedback section
	print $c->h3('Feedback');
	print $CONTACT;
		
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
		}
	}
	print $c->p;
	print $c->submit, $c->defaults('Clear Form');
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

