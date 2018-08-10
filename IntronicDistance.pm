=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

 Stephen Kazakoff <sh.kazakoff@gmail.com>
    
=cut

=head1 NAME

 IntronicDistance

=head1 SYNOPSIS

 mv IntronicDistance.pm ~/.vep/Plugins
 ./vep -i variants.vcf --plugin IntronicDistance

=head1 DESCRIPTION

 A VEP plugin that just returns the HGVS intron start and end offsets.

=cut

package IntronicDistance;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {
    IntronStartOffset => "HGVS intron start offset",
    IntronEndOffset => "HGVS intron end offset",
  };
}

sub run {
  my ($self, $tva) = @_;

  my $hgvs_c = $tva->hgvs_transcript;

  return {} unless $hgvs_c;

  my $start_offset = $tva->hgvs_intron_start_offset;
  my $end_offset = $tva->hgvs_intron_end_offset;

  return {
    "IntronStartOffset" => $start_offset,
    "IntronEndOffset" => $end_offset,
  };
}

1;
