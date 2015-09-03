=head1 LICENSE
                                                                                                                     
 Copyright (c) 1999-2015 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.                                                                      
                                                                                                                     
 This software is distributed under a modified Apache license.                                                       
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html                                                               
                                                                                                                     
=head1 CONTACT                                                                                                       

 Will McLaren <wm2@ebi.ac.uk>
    
=cut

=head1 NAME

 NearestGene

=head1 SYNOPSIS

 mv NearestGene.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --cache --plugin NearestGene

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 finds the nearest gene(s) to a non-genic variant. More than one gene
 may be reported if the genes overlap the variant or if genes are
 equidistant.

 Various parameters can be altered by passing them to the plugin command:

 - limit     : limit the number of genes returned (default: 1)
 - range     : initial search range in bp (default: 1000)
 - max_range : maximum search range in bp (default: 10000)

 Parameters are passed e.g.:

 --plugin NearestGene,limit=3,max_range=50000

=cut

package NearestGene;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my %CONFIG = (
  limit => 1,
  range => 1000,
  max_range => 10000,
);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  
  my $params = $self->params;

  foreach my $param(@$params) {
    my ($key, $val) = split('=', $param);
    die("ERROR: Failed to parse parameter $param\n") unless defined($key) && defined($val);
    die("ERROR: Unknown parameter key $key\n") unless exists($CONFIG{$key});
    $CONFIG{$key} = $val;
  }
  
  return $self;
}

sub feature_types {
  return ['Intergenic','MotifFeature','RegulatoryFeature'];
}

sub variant_feature_types {
  return ['BaseVariationFeature'];
}

sub get_header_info {
  return {
    NearestGene => "Ensembl identifier of nearest gene"
  };
}

sub run {
  my ($self, $vfoa) = @_;
  
  my $vf = $vfoa->base_variation_feature;
  my $loc_string = sprintf("%s:%i-%i", $vf->{chr} || $vf->seq_region_name, $vf->{start}, $vf->{end});
  
  if(!exists($self->{_cache}) || !exists($self->{_cache}->{$loc_string})) {
    $self->{ga} ||= $self->{config}->{ga};
    die("ERROR: Could not get gene adaptor; this plugin does not work in --offline mode\n") unless $self->{ga};
    
    my @result = map {$_->[0]->stable_id} @{
      $self->{ga}->fetch_all_by_outward_search(
        -feature   => $vf,
        -limit     => $CONFIG{limit},
        -range     => $CONFIG{range},
        -max_range => $CONFIG{max_range},
      )
    };
    
    $self->{_cache}->{$loc_string} = scalar @result ? join(",", @result) : undef;
  }
  
  return $self->{_cache}->{$loc_string} ? { NearestGene => $self->{_cache}->{$loc_string} } : {};
}

1;

