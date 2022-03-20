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

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds information from EVE (evolutionary model of variant effect).
 In this plugin, only single nucleotide variants (SNVs) are analyzed,
 multiple nucleotides variants (MNVs) are not included.
 It is only available for GRCh38.

 Please cite EVE publication alongside the VEP if you use this resource:
 https://www.nature.com/articles/s41586-021-04043-8

###################################################
# Bash script to merge all VCFs from EVE dataset. #
###################################################

### BEGIN

# EVE input file can be downloaded from https://evemodel.org/api/proteins/bulk/download/ 
# Input: VCF files by protein (vcf_files_missense_mutations inside zip folder)
# Output: Compressed Merged VCF file (vcf.gz) + index file (.tbi)

DATA_FOLDER='/<PATH-TO>/vcf_files_missense_mutations' # Fill this line
OUTPUT_FOLDER='/<PATH-TO>/eve_plugin' # Fill this line
OUTPUT_NAME='eve_merged.vcf' # Default output name

# Get header from first VCF
cat `ls ${DATA_FOLDER}/*vcf | head -n1` > header

# Get variants from all VCFs and add to a single-file
ls ${DATA_FOLDER}/*vcf | while read VCF; do grep -v '^#' ${VCF} >> variants; done

# Merge Header + Variants in a single file
cat header variants | \
awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > ${OUTPUT_FOLDER}/${OUTPUT_NAME};

# Remove temporary files
rm header variants

# Compress and index
bgzip ${OUTPUT_FOLDER}/${OUTPUT_NAME};

# If not installed, use: sudo apt install tabix
tabix ${OUTPUT_FOLDER}/${OUTPUT_NAME}.gz;

### END

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

  my $assembly = $self->{config}->{assembly};

  die "\nAssembly is not GRCh38, EVE only works with GRCh38. \n" if ($assembly ne "GRCh38");

  my $param_hash = $self->params_to_hash();

  die "\nERROR: No EVE file specified\n" unless defined($param_hash->{file});

  $self->add_file($param_hash->{file});

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;
  return {
    EVE_SCORE => 'Score from EVE model',
    EVE_CLASS   => 'Classification (Benign, Uncertain, or Pathogenic) when setting 75% as uncertain'
  }
}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;

  return {} unless grep {$_->SO_term eq 'missense_variant'} @{$tva->get_all_OverlapConsequences};

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

  # Parsing INFO field
  my ($EVE_SCORE) = $info =~ /EVE=(.*?);/;
  my ($EVE_CLASS) = $info =~ /Class75=(.*?);/;

  return {
    ref => $ref,
    alt => $alt,
    start => $pos,
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

1;
