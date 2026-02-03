=head1 LICENSE
Copyright [2026] EMBL-European Bioinformatics Institute

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
 promotorAI

=head1 SYNOPSIS
 mv promotorAI.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin promotorAI,/path/to/promotorAI.tsv.gz

=head1 DESCRIPTION
 An Ensembl VEP plugin that adds promotorAI scores to promotor variants, predicting their impact on gene expression.

 To download the promotorAI scores file to use with VEP (GRCh38 based),
   Please follow the instructions found in the README at https://github.com/Illumina/PromoterAI.
   You need a valid license agreement as described in the README to obtain and use the promotorAI scores.

 Please cite the promotorAI publication alongside Ensembl VEP if you use this resource:
 https://www.science.org/doi/10.1126/science.ads7373

 Necessary before using the plugin:
   Do the following steps to index the annotations file before using the plugin:
    zcat promotorAI.tsv.gz | sed '1s/.*/#&/' | bgzip > promotorAI.tsv.bgz
    tabix -s 1 -b 2 -e 2 promotorAI.tsv.bgz

 You must have the Bio::DB::HTS module or the tabix utility must be installed
 in your path to use this plugin.
=cut

package promotorAI;

use strict;
use warnings;

use File::Basename;
use File::Spec;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();

  my $file = $self->params->[0];


  my $headers;
  $self->add_file($file);

  open FH, "tabix -fh $file 1:1-1 2>&1 | ";
  while(<FH>){
      next unless /^\#/;
      chomp;
      $_ =~ s/^\#//;
      $self->{headers} = [split];
  }

  close(FH);
  $headers = scalar @{$self->{headers}};

  die("ERROR: Could not find promotorAI scores file\n") unless $file;


  $self->{file_column} = $headers;


  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my %header_info;

  my $vcf_attribute = 'promotorAI_score';
  my $attribute_descr = 'promotorAI score';

  $header_info{ $vcf_attribute } = $attribute_descr;

  return \%header_info;
}

sub run {
  my ($self, $vfoa) = @_;

  my $vf = $vfoa->variation_feature;

  (my $vf_chr = $vf->{chr}) =~ s/^chr//;
  my ($vf_start, $vf_end) = ($vf->{start}, $vf->{end});

  $vf_end = $vf_start if $vf_start > $vf_end;

  my @data = @{ $self->get_data($vf_chr, $vf_start -2, $vf_end) };

  return {} unless @data;

  # Match the relevant alternative allele at the given position
  my $promotorAI_result = {};
  my $vfoa_variation_feature_seq = $vfoa->variation_feature_seq;

  for my $data_candidate (@data) {
    my $candidate_alt_allele = $data_candidate->{'alt'};

    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => [ $vfoa_variation_feature_seq ],
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
        ref  => $data_candidate->{ref},
        alts => [ $candidate_alt_allele ],
        pos  => $data_candidate->{start},
        strand => $data_candidate->{strand}
      }
    );

    if (@$matches){
      $promotorAI_result = $data_candidate->{'result'};
      last;
    }
  }

  return $promotorAI_result;
}

sub parse_data {
  my ($self, $line) = @_;

  my $header = $self->{file_column};
  my ($chr, $pos, $ref, $alt, $gene, $gene_id, $transcript_id, $strand, $tss_position, $score) = split /\t/, $line;

  # VCF-like adjustment of mismatched substitutions for comparison with VEP
  if(length($alt) != length($ref)) {
    my $first_ref = substr($ref, 0, 1);
    my $first_alt = substr($alt, 0, 1);
    if ($first_ref eq $first_alt) {
      $pos++;
      $ref = substr($ref, 1) || "-";
      $alt = substr($alt, 1) || "-";
    }
  }

  my @keys = qw(promotorAI_score);

  my %result;

  @result{@keys} = ($score);

  return {
    ref => $ref,
    alt => $alt,
    start => $pos,
    strand => $strand,
    gene_id => $gene_id,
    transcript_id => $transcript_id,
    result => \%result
  }
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;

