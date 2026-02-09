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
 ./vep -i variations.vcf --plugin gnomADMt,file=/path/to/gnomad-chrM.vcf.gz

=head1 DESCRIPTION
 An Ensembl VEP plugin that retrieves allele frequencies for Mitochnodrial variants
 from the gnomAD genomes mitochnodrial annotations files, available here:
   https://gnomad.broadinstitute.org/downloads#v3-mitochondrial-dna

 Options are passed to the plugin as key=value pairs:
   file : (mandatory) Tabix-indexed VCF file from gnomAD
   fields : Colon-separated list of information from mitochnodrial variants to
            output (default: 'AC_hom:AC_het:AF_hom:AF_het:AN:max_hl');
            keyword 'all' can be used to print all fields;
            Available fields are the names of INFO fields.

 To download the gnomad Mitochondrial allele annotations file in VCF format (GRCh38 based) and the corresponding index file:
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz --no-check-certificate
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz.tbi --no-check-certificate

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

sub _get_valid_fields {
  my $selected  = shift;
  my $available = shift;

  # return all available fields when using 'all'
  return $available if $selected eq 'all';
  my @fields = split(/:/, $selected);

  # check if the selected fields exist
  my @valid;
  my @invalid;
  for my $field (@fields) {
    if ( grep { $_ eq $field } @$available ) {
      push(@valid, $field);
    } else {
      push(@invalid, $field);
    }
  }

  die "ERROR: all fields given are invalid. Available fields are:\n" .
    join(", ", @$available)."\n" unless @valid;
  warn "gnomADMt plugin: WARNING: the following fields are not valid and were ignored: ",
    join(", ", @invalid), "\n" if @invalid;

  return \@valid;
}

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();
  my $params = $self->params_to_hash();

  my $file = $params->{file};

  # get INFO fields names from VCF
  my $vcf_file = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($file);
  my $info = $vcf_file->get_metadata_by_pragma('INFO');
  my $info_ids = [ map { $_->{ID} } @$info ];
  my $info_descriptions = { map { $_->{ID} => $_->{Description} } @$info };

  # check if AF_hom and AF_het fields exists
  if (!defined $info_descriptions->{AF_hom} || !defined $info_descriptions->{AF_het}) {
    die "ERROR: Provided file does not contain AF_hom or AF_het fields. Please provide the file as downloaded from GnomAD.\n"
  }

  $self->add_file($file);

  # Process the requested fields
  my $requested_fields_str = defined $params->{fields} ? $params->{fields} : 'AC_hom:AC_het:AF_hom:AF_het:AN:max_hl';
  my $fields = _get_valid_fields($requested_fields_str, $info_ids);

  $self->{fields} = $fields;
  $self->{fields_descriptions} = { map { $_ => $info_descriptions->{$_} } @$fields };

  my $prefix = 'gnomAD';
  $self->{prefix} = $prefix;

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my $prefix = $self->{prefix} . '_';
  my %header_info = map {
    ''.$prefix.$_ => $self->{fields_descriptions}->{$_}
  } @{$self->{fields}};

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

sub _parse_vcf {
  my ($self, $line) = @_;

  my ($chrom, $start, $id, $ref, $alt, $qual, $filter, $info) = split /\t/, $line;
  # VCF-like adjustment of mismatched substitutions for comparison with VEP
  if(length($alt) != length($ref)) {
    my $first_ref = substr($ref, 0, 1);
    my $first_alt = substr($alt, 0, 1);
    if ($first_ref eq $first_alt) {
      $start++;
      $ref = substr($ref, 1) || "-";
      $alt = substr($alt, 1) || "-";
    }
  }

  my %data = (
    'chromosome' => $chrom,
    'start'      => $start,
    'identifier' => $id,
    'alleles'    => $ref.'/'.$alt,
    'ref'        => $ref,
    'alt'        => $alt,
    'quality'    => $qual,
    'filter'     => $filter,
  );
  # fetch data from all INFO fields
  for my $field ( split /;/, $info ) {
    my ($key, $value) = split /=/, $field;
    $data{$key} = $value;
  }

  return \%data;
}

sub parse_data {
  my ($self, $line) = @_;

  my $vcf_data = $self->_parse_vcf($line);

  my %keys = map {
    $_ => join('_', $self->{prefix}, $_)
  }
  @{$self->{fields}};

  # Filter VCF data down to selected fields
  # and prefix the field names
  my $result = {map { $keys{$_} => $vcf_data->{$_} } @{$self->{fields}}};

  return {
    ref => $vcf_data->{ref},
    alt => $vcf_data->{alt},
    start => $vcf_data->{start},
    result => $result
  }
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;

