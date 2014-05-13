=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and                                                   
 Genome Research Limited.  All rights reserved.                                                                      

 This software is distributed under a modified Apache license.                                                       
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html                                                               

=head1 CONTACT                                                                                                       

 Will McLaren <wm2@ebi.ac.uk>
    
=cut

=head1 NAME

 CADD

=head1 SYNOPSIS

 mv miRNA.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin miRNA

=head1 DESCRIPTION

 A VEP plugin that determines where in the secondary structure of a miRNA a
 variant falls.
 
=cut

package miRNA;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;
  return {
    miRNA => 'SO term for miRNA component containing the variant'
  }
}

sub run {
  my ($self, $tva) = @_;
  
  my $tv = $tva->transcript_variation;
  my $tr = $tva->transcript;
  
  # obviously this only works for *RNA transcripts
  return {} unless $tr->biotype =~ /RNA/;
  
  # get attribute if already cached
  my ($attrib) = @{$tr->get_all_Attributes('ncRNA')};
  
  # bit of a cheat to get attrib if ncRNA attribute hasn't been cached
  if(!$attrib && defined($self->{config}->{ta})) {
    delete $tr->{attributes};
    $tr->{adaptor} = $self->{config}->{ta};
    ($attrib) = @{$tr->get_all_Attributes('ncRNA')};
  }
  
  return {} unless $attrib;
  
  # split out string to get coords and structure string
  my ($start, $end, $struct) = split /\s+|\:/, $attrib->value;
  return {} unless $struct && $struct =~ /[\(\.\)]+/;
  
  # variant not in given structure?
  return { miRNA => 'None' } unless $tv->cdna_start <= $end && $tv->cdna_end >= $start;
  
  # parse out structure
  my @struct;
  while($struct =~ m/([\.\(\)])([0-9]+)?/g) {
    my $num = $2 || 1;
    push @struct, $1 for(1..$num);
  }
  
  # get struct element types overlapped by variant
  my %chars;
  for my $pos($tv->cdna_start..$tv->cdna_end) {
    $pos -= $start;
    next if $pos < 0 or $pos > scalar @struct;
    $chars{$struct[$pos]} = 1;
  }
  
  # map element types to SO terms
  my %map = (
    '(' => 'miRNA_stem',
    ')' => 'miRNA_stem',
    '.' => 'miRNA_loop'
  );
  
  return {
    miRNA => join(",", sort map {$map{$_}} keys %chars)
  };
}

1;

