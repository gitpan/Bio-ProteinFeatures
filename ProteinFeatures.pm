# $Id: ProteinFeatures.pm,v 1.5 2003/09/08 09:26:14 cvspub Exp $
package Bio::ProteinFeatures;
use strict;
our $VERSION = '0.01';

=pod

=head1 NAME

Bio::ProteinFeatures - Deriving features of amino acid sequences

=head1 SYNOPSIS

  use Bio::ProteinFeatures;

  $pf = new Bio::ProteinFeatures;

  $pf->sequence($sequence_string);


  # you may use Data::Dumper to see the result.
  use Data::Dumper;
  print Dumper $pf->features();


=head1 DESCRIPTION

This module applies several statistical methods on amino acid sequences for deriving various useful features for identifying sequences and they may be used to measure similarities between sequences. You may also use this module to do coarse matching before doing Blast.

=head1 METHODS

=cut


# ----------------------------------------------------------------------

=pod

=head2 new

You can set the sequence on invoking the constructor.

    $pf = new Bio::ProteinFeatures(sequence => $sequence_string);

Or set it using the next method.

=cut

sub new {
    my $pkg = shift;
    my %arg = @_;
    bless {
	_seq => uc $arg{sequence},
	_seqlen => length $arg{sequence},
    }, 
}


# ----------------------------------------------------------------------

=pod

=head2 sequence

Set or get the sequence string

    # set the sequence
    $pf->sequence($sequence);

    # return the sequence
    $pf->sequence();

=cut

sub sequence {
    my $this = shift;
    if( $_[0] ){
	$this->{_seq} = uc shift();
	$this->{_seqlen} = length $this->{_seq};
    }
    $this->{_seq};
}


# ----------------------------------------------------------------------

=pod

=head2 features

The features this module deals with are listed below.

=cut

use constant POLAR => 'polar';
use constant NEUTRAL => 'neutral';
use constant HYDROPHOBE => 'hydrophobic';
our @aa = qw(A C D E F G H I K L M N P Q R S T V W Y);

our %category = ( R => POLAR, K => POLAR, E => POLAR, D => POLAR, Q =>
		 POLAR, N => POLAR, G => NEUTRAL, A => NEUTRAL, S =>
		 NEUTRAL, T => NEUTRAL, P => NEUTRAL, H => NEUTRAL, Y
		 => NEUTRAL, C => HYDROPHOBE, V => HYDROPHOBE, L =>
		 HYDROPHOBE, I => HYDROPHOBE, M => HYDROPHOBE, F =>
		 HYDROPHOBE, W => HYDROPHOBE
		 );

use List::Util qw(sum reduce);
sub features {
    my $this = shift;
    my %features;

=pod

=head3 composition

Amino acids are grouped into three categories: polar, neutral, and hydrophobic. The methods calculates the compositions of the three groups of amino acids.

=cut

    my %composition;
    for( my $i ; $i<$this->{_seqlen} ; $i++){
	$composition{$category{substr($this->{_seq}, $i, 1)}}++;
    }
    %composition = map {$_ => $composition{$_}/$this->{_seqlen}} keys %composition;
    $features{composition} = \%composition;


=pod

=head3 transition probability

Characterizes the percent frequency with which group A is followed by group B or B is followed by A.

=cut

    my %t_freq;
    for( my $i ; $i<$this->{_seqlen}-1 ; $i++){
	next if $category{substr($this->{_seq}, $i, 1)} eq
	    $category{substr($this->{_seq}, $i+1, 1)};
	$t_freq{join q/_/, sort $category{substr($this->{_seq}, $i, 1)},
		    $category{substr($this->{_seq}, $i+1, 1)}}++;
    }
    $features{transition_probability} = {
	map{$_=>$t_freq{$_}/($this->{_seqlen}-1)} keys %t_freq };

=pod

=head3 accumulative distribution

Sequences are cut into 5 sections. It calculates the accumulative probabilities of a certain group within a section.

=cut

    my %dist;
    foreach my $sec (1..5){
	my $tmpend = $this->{_seqlen}*($sec/5);
	$tmpend = $tmpend-1 if $tmpend == $this->{_seqlen};
	foreach (0..$tmpend){
	    $dist{$category{substr($this->{_seq}, $_, 1)}}->{$sec}++;
	}
    }
    foreach my $cat (keys %dist){
	foreach my $sec (1..5){
	    $dist{$cat}->{$sec} /= $dist{$cat}->{5};
	}
    }
    $features{accumulative_distribution} = \%dist;

=pod

=head3 per-amino-acid probability    

Calculates per-se probability of each amino acid.

=cut

    my %aafreq;
    for( my $i ; $i<$this->{_seqlen} ; $i++){
	$aafreq{substr($this->{_seq}, $i, 1)}++;
    }
    my %aaprob = map{$_=>$aafreq{$_}/$this->{_seqlen}} keys %aafreq;
    $features{aa_prob} = \%aaprob;

=pod

=head3 first order energy

summation prob(i**2) for each i of amino acids.

=cut

    my $first_order_energy = 
	reduce { $a + $b } map{$_**2} values %aaprob;
    $features{first_order_energy} = $first_order_energy;

=pod

=head3 first order entropy

summation -prob(i)*log(prob(i)) for each i of amino acids.

=cut

    my $first_order_entropy = 
	reduce { $a + $b } map{-$_*log$_} values %aaprob;
    $features{first_order_entropy} = $first_order_entropy;


=pod

=head3 histogram difference

Calculates the difference of the numbers of two neighboring amino acids.


=cut

    my ($histo_diff);

    foreach (my $i=0; $i<@aa-1; $i++){
	my $diff = abs($aaprob{$aa[$i]}-$aaprob{$aa[$i+1]});
	### histogram difference
	$histo_diff += $diff;
	
    }
    $features{histogram_difference} = $histo_diff;

=pod

=head3 AA pair probability

Probabilities of amino acid bigrams.

=cut

    my (%pair_prob);
    foreach (my $i=0; $i<@aa-1; $i++){
	$pair_prob{substr($this->{_seq}, $i, 2)}++;
    }
    %pair_prob = map{$_=>$pair_prob{$_}/($this->{_seqlen}-1)} keys %pair_prob;
    my $second_order_energy = reduce{$a+$b} map{$_**2} values %pair_prob;
    my $second_order_entropy = reduce{$a+$b} map{-$_*log$_} values %pair_prob;

=pod

=head3 average seperation between two amino acid of the same group

Counts the average number of characters between two amino acids of the same group.

=cut

    my $avg_sep = 0;
    foreach my $aa (@aa){
	my $left = 0;
	my $right = 0;
	$left = index($this->{_seq}, $aa);
	my $local_sep_sum = 0;
	while($left!=-1){
	    if(($right = index($this->{_seq}, $aa, $left+1)) != -1){
		$local_sep_sum += $right-$left;
	    }
	    $left = $right;
	}
	$avg_sep += $aafreq{$aa} ?
	    $local_sep_sum/($this->{_seqlen}*$aafreq{$aa}) : 0;
    }

    \%features;
}





1;
__END__

=head1 COPYRIGHT

xern E<lt>xern@cpan.orgE<gt>

This module is free software; you can redistribute it or modify it under the same terms as Perl itself.

=cut
