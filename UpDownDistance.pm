=head1 LICENSE
                                                                                                                     
 Copyright (c) 1999-2013 The European Bioinformatics Institute and                                                   
 Genome Research Limited.  All rights reserved.                                                                      
                                                                                                                     
 This software is distributed under a modified Apache license.                                                       
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html                                                               
                                                                                                                     
=head1 CONTACT                                                                                                       

 William McLaren <wm2@ebi.ac.uk>
    
=cut

=head1 NAME

 UpDownDistance

=head1 SYNOPSIS

 mv UpDownDistance.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin UpDownDistance,10000

=head1 DESCRIPTION

 Change the distance used by the VEP to call upstream and downstream consequence
 types. Defaults to 5000. Up- and downstream distances can be set separately by
 passing two arguments:
 
 --plugin UpDownDistance,[up],[down]
=cut

package UpDownDistance;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);
  
  # change up/down
  my $up = $self->params->[0] || 5000;
  my $down = $self->params->[1] || $up;
  $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE = $up;
  $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE = $down;
  
  return $self;
}

sub run { return {} };

1;

