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

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

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
  my $alt_allele = $tva->variation_feature_seq;
  my $ref_allele = $vf->ref_allele_string;

  return {} unless $alt_allele =~ /^[ACGT-]+$/;

  my @data = @{
    $self->get_data(
      $vf->{chr},
      $vf->{start} - 2,
      $vf->{end}
    )
  };

  return {} unless(@data);

  foreach my $variant (@data) {
    # This check should disappear with multiple snps codons
    return {} unless($variant->{start} and $variant->{ref} and $variant->{alt});

    my $matches = get_matched_variant_alleles(
      {
        ref    => $ref_allele,
        alts   => [$alt_allele],
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
       ref  => $variant->{ref},
       alts => [$variant->{alt}],
       pos  => $variant->{start},
      }
    );

    return $variant->{result} if (@$matches);
  }

  return {};

}

sub parse_data {
  my ($self, $line) = @_;

  # Parsing VCF fields
  my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info) = split /\t/, $line; 

  my $allele_difference = ($ref ^ $alt) =~ tr/\0//c;

  # TODO Prepare multiple SNPs codons part
  if ($allele_difference > 1){
    return {};
  }

  # Parsing INFO field
  my ($EVE_SCORE) = $info =~ /EVE=(.*?);/;
  my ($EVE_CLASS) = $info =~ /Class70=(.*?);/;

  my %new_snp = _get_snp_from_codon($ref, $alt, $pos);

  return {
    ref => $new_snp{ref},
    alt => $new_snp{alt},
    start => $new_snp{pos},
    result => {
      EVE_SCORE   => $EVE_SCORE,
      EVE_CLASS   => $EVE_CLASS
    }
  }
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

# Additional function
sub _get_snp_from_codon {
  my ($ref, $alt, $pos) = @_;

  my $mask = $ref ^ $alt;
  while ($mask =~ /[^\0]/g) {
      $ref = substr($ref,$-[0],1);
      $alt = substr($alt,$-[0],1);
      $pos += $-[0];
  }

  return (
    ref => $ref,
    alt => $alt,
    pos => $pos
  )

}

1;
