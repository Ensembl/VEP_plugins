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

 CAPICE

=head1 SYNOPSIS

 mv CAPICE.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin CAPICE,/FULL_PATH_TO_CAPICE_FILE/capice_v1.0_build37_snvs.tsv.gz,/FULL_PATH_TO_CAPICE_FILE/capice_v1.0_build37_indels.tsv.gz
 ./filter_vep -i variant_effect_output.txt --filter "CAPICE_SCORE >= 0.02"

=head1 DESCRIPTION

 A VEP plugin that retrieves CAPICE scores for variants from one or more
 tabix-indexed CAPICE data files, in order to predict their pathogenicity.
 
 Please cite the CAPICE publication alongside the VEP if you use this resource:
 https://pubmed.ncbi.nlm.nih.gov/32831124/
 
 The tabix utility must be installed in your path to use this plugin. The CAPICE
 data files can be downloaded from https://zenodo.org/record/3928295

 The plugin works with CAPICE files for GRCh37. The plugin only reports scores
 and does not consider any additional annotations from a CAPICE file. It is
 therefore sufficient to use CAPICE files without the additional annotations.

 To filter results, please use filter_vep with the output file or standard
 output. Documentation on filter_vep is available at:
 https://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html

 For recommendations on threshold selection, please read the CAPICE publication.

=cut

package CAPICE;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();
  
  # Check files in arguments
  my @params = @{$self->params};
  die "ERROR: No CAPICE files specified" unless @params > 0;
  $self->add_file($_) for @params;

  return $self;
}

sub feature_types {
  return ['Feature','Intergenic'];
}

sub get_header_info {
  return {
    CAPICE_SCORE => 'CAPICE score to predict the pathogenicity of SNVs and InDels'
  }
}

sub run {
  my ($self, $tva) = @_;
  
  my $vf = $tva->variation_feature;
  
  # get allele
  my $allele = $tva->variation_feature_seq;
  
  return {} unless $allele =~ /^[ACGT-]+$/;

  my @data = @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})};

  foreach (@data) {
    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => [$allele],
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
       ref  => $_->{ref},
       alts => [$_->{alt}],
       pos  => $_->{start},
      }
    );
    return $_->{result} if (@$matches);
  }
  return {};
}

sub parse_data {
  my ($self, $line) = @_;
  my ($chrom, $start, $ref, $alt, $score) = split /\t/, $line;

  # VCF-like adjustment of mismatched substitutions for comparison with VEP
  # my $end = ($start + length($ref)) - 1;
  if(length($alt) != length($ref)) {
    my $first_ref = substr($ref, 0, 1);
    my $first_alt = substr($alt, 0, 1);
    if ($first_ref eq $first_alt) {
      $start++;
      $ref = substr($ref, 1);
      $alt = substr($alt, 1);
      $ref ||= '-';
      $alt ||= '-';
    }
  }
  return {
    ref => $ref,
    alt => $alt,
    start => $start,
    # end => $end,
    result => {
      CAPICE_SCORE => $score
    }
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
