=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

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

 Paralogues

=head1 SYNOPSIS

 mv Paralogues.pm ~/.vep/Plugins

 # Fetch Ensembl variants in paralogue proteins (requires database access) -- slowest mode
 ./vep -i variations.vcf --database --plugin Paralogues

 # Fetch variants from custom VCF file based on Ensembl paralogues (requires database access)
 ./vep -i variations.vcf --database --plugin Paralogues,vcf=/path/to/file.vcf

=head1 DESCRIPTION

 A VEP plugin that fetches variants from paralogue proteins.
 
 The tabix utility must be installed in your path to use this plugin.

=cut

package Paralogues;

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
  my $params = $self->params_to_hash();
  $self->add_file($params->{vcf}) if defined $params->{vcf};

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my $header = {
    PARALOGUE_VARIANT_ID          => 'Paralogue variant identifier',
    PARALOGUE_VARIANT_SEQNAME     => 'Paralogue variant sequence region name',
    PARALOGUE_VARIANT_START       => 'Paralogue variant start',
    PARALOGUE_VARIANT_ALLELES     => 'Paralogue variant alleles',
  };

  if ( @{$self->{_files}} ) {
    $header->{'PARALOGUE_VARIANT_INFO_*'} = 'Paralogue variant fields based on VCF INFO field';
  } else {
    $header->{'PARALOGUE_VARIANT_STRAND'}      = 'Paralogue variant strand';
    $header->{'PARALOGUE_VARIANT_END'}         = 'Paralogue variant end';
    $header->{'PARALOGUE_VARIANT_CLNSIG'}      = 'Paralogue variant clinical significance';
    $header->{'PARALOGUE_VARIANT_SOURCE'}      = 'Paralogue variant source';
    $header->{'PARALOGUE_VARIANT_CONSEQUENCE'} = 'Paralogue variant consequence type(s)';
  }
  return $header;
}

sub _join_results {
  my $self = shift;
  my $all_results = shift;
  my $res = shift;

  if ($self->{config}->{output_format} eq 'json' || $self->{config}->{rest}) {
    # Group results for JSON and REST
    $all_results->{"Paralogues"} = [] unless defined $all_results->{"Paralogues"};
    push(@{ $all_results->{"Paralogues"} }, $res);
  } else {
    # Create array of results per key
    for (keys %$res) {
      $all_results->{$_} = [] if !$all_results->{$_};
      push(@{ $all_results->{$_} }, $res->{$_} || "NA");
    }
  }
  return $all_results;
}

sub _get_transcript_homologies {
  my $self       = shift;
  my $transcript = shift;

  my $config  = $self->{config};
  my $species = $config->{species};
  my $reg     = $config->{reg};

  # reconnect to DB without species param
  if ($config->{host}) {
    $reg->load_registry_from_db(
      -host       => $config->{host},
      -user       => $config->{user},
      -pass       => $config->{password},
      -port       => $config->{port},
      -db_version => $config->{db_version},
      -no_cache   => $config->{no_slice_cache},
    );
  }
  my $ha = $reg->get_adaptor( "multi", "compara", "homology" );

  $self->{ga} ||= $reg->get_adaptor($species, 'core', 'gene');
  my $gene = $self->{ga}->fetch_by_stable_id( $transcript->{_gene_stable_id} );
  my $homologies = $ha->fetch_all_by_Gene(
    $gene, -METHOD_LINK_TYPE => 'ENSEMBL_PARALOGUES', -TARGET_SPECIES => $species);
  return $homologies;
}

sub _get_paralogue_coords {
  my $self     = shift;
  my $tva      = shift;
  my $homology = shift;

  my $aln = $homology->get_SimpleAlign;
  my $translation_id    = $tva->transcript->translation->stable_id;
  my $translation_start = $tva->base_variation_feature_overlap->translation_start;

  # identify paralogue protein
  my @proteins = keys %{ $aln->{'_start_end_lists'} };
  my $paralogue;
  if ($translation_id eq $proteins[0]) {
    $paralogue = $proteins[1];
  } elsif ($translation_id eq $proteins[1]) {
    $paralogue = $proteins[0];
  } else {
    return;
  }

  # get genomic coordinates from paralogue
  my $col    = $aln->column_from_residue_number($translation_id, $translation_start);
  my $coords = $aln->get_seq_by_id($paralogue)->location_from_column($col);
  return unless defined $coords and $coords->location_type eq 'EXACT';

  my $config  = $self->{config};
  my $species = $config->{species};
  my $reg     = $config->{reg};

  $self->{ta}     ||= $reg->get_adaptor($species, 'core', 'translation');
  my $para_tr       = $self->{ta}->fetch_by_stable_id($paralogue)->transcript;
  my ($para_coords) = $para_tr->pep2genomic($coords->start, $coords->start);

  my $chr   = $para_tr->seq_region_name;
  my $start = $para_coords->start;
  my $end   = $para_coords->end;
  return ($chr, $start, $end);
}

sub _summarise_vf {
  my $vf = shift;
  return {
    PARALOGUE_VARIANT_ID          => $vf->name,
    PARALOGUE_VARIANT_SEQNAME     => $vf->seq_region_name,
    PARALOGUE_VARIANT_STRAND      => $vf->strand,
    PARALOGUE_VARIANT_START       => $vf->seq_region_start,
    PARALOGUE_VARIANT_END         => $vf->seq_region_end,
    PARALOGUE_VARIANT_ALLELES     => $vf->allele_string,
    PARALOGUE_VARIANT_CLNSIG      => join('/', @{$vf->get_all_clinical_significance_states}),
    PARALOGUE_VARIANT_SOURCE      => $vf->source_name,
    PARALOGUE_VARIANT_CONSEQUENCE => $vf->display_consequence,
  };
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  my $transcript = $tva->transcript;
  my $translation_start = $tva->base_variation_feature_overlap->translation_start;
  return {} unless defined $translation_start;

  my $config  = $self->{config};
  my $reg     = $config->{reg};
  my $species = $config->{species};

  my $all_results = {};
  my $homologies = $self->_get_transcript_homologies($transcript);
  for my $homology (@$homologies) {
    my ($chr, $start, $end) = $self->_get_paralogue_coords($tva, $homology);
    next unless defined $chr and defined $start and defined $end;

    if (@{$self->{_files}}) {
      # get variants from custom VCF file
      my @data = @{$self->get_data($chr, $start - 2, $end)};
      $all_results = $self->_join_results($all_results, $_) for @data;
    } else {
      # get Ensembl variants from mapped genomic coordinates
      $self->{sa}  ||= $reg->get_adaptor($species, 'core', 'slice');
      $self->{vfa} ||= $reg->get_adaptor($species, 'variation', 'variationfeature');

      my $slice = $self->{sa}->fetch_by_region('chromosome', $chr, $start, $end);
      foreach my $var ( @{ $self->{vfa}->fetch_all_by_Slice($slice) } ) {
        my $res = _summarise_vf($var);
        $all_results = $self->_join_results($all_results, $res);
      }
    }
  }
  return $all_results;
}

sub parse_data {
  my ($self, $line) = @_;
  my ($chrom, $start, $id, $ref, $alt, $qual, $filter, $info) = split /\t/, $line;

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

  my $data = {
    PARALOGUE_VARIANT_SEQNAME => $chrom,
    PARALOGUE_VARIANT_START   => $start,
    PARALOGUE_VARIANT_ID      => $id,
    PARALOGUE_VARIANT_ALLELES => "$ref/$alt"
  };

  # include data from INFO fields
  for my $field ( split /;/, $info ) {
    my ($key, $value) = split /=/, $field;
    $data->{'PARALOGUE_VARIANT_INFO_' . $key} = $value;
  };
  return $data;
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
