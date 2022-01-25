=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>

=cut

=head1 NAME

 EVE

=head1 SYNOPSIS

 cp EVE.pm ${HOME}/.vep/Plugins
 ./vep -i variations.vcf --plugin EVE,file=/path/to/eve/data.vcf.gz

=head1 DESCRIPTION



=cut
package EVE;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  my $param_hash = $self->params_to_hash();

  $self->add_file($param_hash->{file});

  return $self;
}

sub feature_types {
  return ['Feature'];
}

sub get_header_info {
  my $self = shift;
  return {
    EVE_SCORE => 'Score from EVE model',
    EVE_CLASS   => 'Classification (Benign, Uncertain, or Pathogenic) when setting 70% as uncertain'
  }
}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;

  # get allele
  my $allele = $tva->variation_feature_seq;

  return {} unless $allele =~ /^[ACGT-]+$/;

  my @data = @{
    $self->get_data(
      $vf->{chr},
      $vf->{start} - 2,
      $vf->{end}
    )
  };

  return {} unless(@data);

  foreach my $variant (@data) {
    my @result;
    
    my $EVE_SCORE = $variant->{EVE_SCORE};
    my $EVE_CLASS = $variant->{EVE_CLASS};

    return {
      EVE_SCORE   => $EVE_SCORE,
      EVE_CLASS => $EVE_CLASS
    };

  }

}

sub parse_data {
  my ($self, $line) = @_;

  # Data is a VCF file:
  # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
  # 1       961387  .       GCC     AAA     .       .       EVE=0.08912364373354463;EnsTranscript=ENST0....
  # 1       961387  .       GCC     AAC     .       .       EVE=0.08739485410893585;EnsTranscript=ENST0

  # Parsing VCF fields
  my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info) = split /\t/, $line; 

  my $allele_difference = ($ref ^ $alt) =~ tr/\0//c;

  if ($allele_difference > 1){
    return {};
  }

  # Parsing INFO field
  my ($EVE_SCORE) = $info =~ /EVE=(.*?);/;
  my ($EVE_CLASS) = $info =~ /Class70=(.*?);/;

  return {
    chrom => $chrom,
    ref => $ref,
    alt => $alt,
    start => $pos,
    EVE_SCORE   => $EVE_SCORE,
    EVE_CLASS   => $EVE_CLASS
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
