=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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
 # Download CAPICE SNVs, InDels and index (TBI) files to the same path
 # - capice_v1.0_build37_indels.tsv.gz
 # - capice_v1.0_build37_indels.tsv.gz.tbi
 # - capice_v1.0_build37_snvs.tsv.gz
 # - capice_v1.0_build37_snvs.tsv.gz.tbi
 ./vep -i variations.vcf --plugin CAPICE,snv=/FULL_PATH_TO_CAPICE_FILE/capice_v1.0_build37_snvs.tsv.gz,indels=/FULL_PATH_TO_CAPICE_FILE/capice_v1.0_build37_indels.tsv.gz
 ./filter_vep -i variant_effect_output.txt --filter "CAPICE_SCORE >= 0.02"

=head1 DESCRIPTION

 A VEP plugin that retrieves CAPICE scores for variants from one or more
 tabix-indexed CAPICE data files, in order to predict their pathogenicity.
 
 Please cite the CAPICE publication alongside the VEP if you use this resource:
 https://pubmed.ncbi.nlm.nih.gov/32831124/
 
 The tabix utility must be installed in your path to use this plugin. The CAPICE
 SNVs, InDels and respective index (TBI) files for GRCh37 can be downloaded from
 https://zenodo.org/record/3928295

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
  my @files; 
  my $params = $self->params_to_hash();

  if (!keys %$params){
    @files = @{$self->params};
  } else {
    for my $key (keys %{$params}) {
      push @files, $params->{$key};
    }
  }

  die "ERROR: No CAPICE files specified\n" unless @files > 0;
  $self->add_file($_) for @files;

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
  my $allele = $tva->base_variation_feature->alt_alleles;

  my @data = @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})};

  foreach (@data) {
    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => $allele,
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
