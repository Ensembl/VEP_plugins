=head1 LICENSE
Copyright [2025-2026] EMBL-European Bioinformatics Institute

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
 gnomADMt

=head1 SYNOPSIS
 mv gnomADMt.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin gnomADMt,/path/to/gnomad-chrM-reduced-annotations.tsv.gz

=head1 DESCRIPTION
 An Ensembl VEP plugin that retrieves allele frequencies for Mitochnodrial variants
 from the gnomAD genomes mitochnodrial (reduced) annotations files, available here:
   https://gnomad.broadinstitute.org/downloads

 To download the gnomad Mitochondrial allele annotations file in TSV format (GRCh38 based):
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv --no-check-certificate

 Necessary before using the plugin
   The following steps are necessary to tabix the gnomad annotations file :
    sed '1s/.*/#&/' gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv | bgzip > gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv.bgz
    tabix -s 1 -b 2 -e 2 gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv.bgz

 If you use this plugin, please see the terms and data information:
   https://gnomad.broadinstitute.org/terms

 You must have the Bio::DB::HTS module or the tabix utility must be installed
 in your path to use this plugin.
=cut

package gnomADMt;

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


  my $prefix = 'gnomAD';
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

  die("ERROR: Could not find any $prefix coverage files\n") unless $file;


  $self->{prefix} = $prefix;
  $self->{file_column} = $headers;


  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my $prefix = $self->{prefix};
  my %header_info;


  for (qw(AC_hom AC_het AF_hom AF_het AN max_observed_heteroplasmy)) {
    my $vcf_attribute = join('_', $prefix, $_);
    my $attribute_descr = '';

    if ( $_ =~ /^AF_/ ) {
      $attribute_descr .= 'allele frequency';
    }
    elsif ( $_ =~ /^AC_/ ) {
      $attribute_descr .= 'allele count';
    }

    if ( $_ =~ /_hom$/ ) {
      $attribute_descr .= 'in homozygous carriers';
    }
    elsif ( $_ =~ /_het$/ ) {
      $attribute_descr .= 'in heterozygous carriers';
    }

    if ( $_ eq 'AN' ) {
      $attribute_descr = 'no. of callable samples at this site';
    }

    $header_info{ $vcf_attribute } = $attribute_descr;
  }

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
  my $gnomad_freqs = {};
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
      }
    );

    if (@$matches){
      $gnomad_freqs = $data_candidate->{'result'};
      last;
    }
  }

  return $gnomad_freqs;
}

sub parse_data {
  my ($self, $line) = @_;

  my $prefix = $self->{prefix};
  my $header = $self->{file_column};
  my ($chr, $pos, $ref, $alt, $filters, @af) = split /\t/, $line;

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

  my @keys;

  @keys = map {
    join('_', $prefix, $_)
  }
  qw(filters AC_hom AC_het AF_hom AF_het AN max_observed_heteroplasmy);

  my %result;

  @result{@keys} = ($filters, @af);

  return {
    ref => $ref,
    alt => $alt,
    start => $pos,
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

