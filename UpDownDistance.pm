=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

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

