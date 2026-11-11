=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2026] EMBL-European Bioinformatics Institute

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

AlphaGenomeSplicing - annotate variants with AlphaGenome merged splicing scores

=head1 SYNOPSIS

 mv AlphaGenomeSplicing.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin AlphaGenomeSplicing,file=/FULL_PATH/combined_alphagenome_splicing_snvs.tsv.gz

=head1 DESCRIPTION

 An Ensembl VEP plugin that retrieves AlphaGenome merged splicing scores for single nucleotide variants
 from a tabix-indexed, bgzip-compressed TSV file.

 The scores were developed by Google DeepMind and combine AlphaGenome's splice site, splice site usage and
 splice junction predictions into one score per variant. Higher scores mean a larger predicted effect on
 splicing. They are available for GRCh38 SNVs only, for non-commercial use only.

 The data is "AlphaGenome SNV merged splicing scores" from :
 https://deepmind.google.com/science/alphagenome/downloads

 Download combined_splicing_snvs_tabix.zip and extract combined_alphagenome_splicing_snvs.tsv.gz and its
 .tbi index, which are already bgzip-compressed and tabix-indexed. The file has the columns:

 #CHROM POS REF ALT alphagenome_splicing

 The score is reported as AlphaGenome_splicing.

=cut

package AlphaGenomeSplicing;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my $SCORE_NAME = 'AlphaGenome_splicing';
my $SCORE_DESCRIPTION = 'AlphaGenome merged splicing score (Google DeepMind), combining splice site, splice site usage and splice junction predictions. Higher => larger predicted effect on splicing. Non-commercial use only.';

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();

  my $params = $self->params_to_hash();
  my ($positional_file) = grep { index($_, '=') == -1 } @{$self->params};
  my $file = $params->{file} || $positional_file;

  die "ERROR: No AlphaGenome splicing file specified\n" unless $file;

  $self->add_file($file);

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  return { $SCORE_NAME => $SCORE_DESCRIPTION };
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  my ($start, $end) = ($vf->{start}, $vf->{end});
  my $ref = $vf->ref_allele_string;
  my $allele = $tva->variation_feature_seq;

  return {} unless defined($start) && defined($end) && $start == $end;
  return {} unless defined($ref) && $ref =~ /^[ACGT]$/;
  return {} unless defined($allele) && $allele =~ /^[ACGT]$/;

  foreach my $record (@{$self->get_data($vf->{chr}, $start, $end)}) {
    my $matches = get_matched_variant_alleles(
      {
        ref    => $ref,
        alts   => [$allele],
        pos    => $start,
        strand => $vf->strand,
      },
      {
        ref  => $record->{ref},
        alts => [$record->{alt}],
        pos  => $record->{start},
      }
    );

    return { $SCORE_NAME => $record->{score} } if @{$matches};
  }

  return {};
}

sub parse_data {
  my ($self, $line) = @_;
  my ($chr, $start, $ref, $alt, $score) = split /\t/, $line, 6;

  return unless defined($score);
  chomp $score;

  return {
    chr   => $chr,
    start => $start,
    end   => $start,
    ref   => $ref,
    alt   => $alt,
    score => $score,
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
