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

 DosageSensitivity

=head1 SYNOPSIS

 mv DosageSensitivity.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin DosageSensitivity,file=/FULL_PATH_TO/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz
 ./vep -i variations.vcf --plugin DosageSensitivity,file=/FULL_PATH_TO/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz,cover=1

=head1 DESCRIPTION

 A VEP plugin that retrieves haploinsufficiency and triplosensitivity probability scores
 for affected genes from a dosage sensitivity catalogue published in paper -
 https://www.sciencedirect.com/science/article/pii/S0092867422007887
 
 Please cite the above publication alongside the VEP if you use this resource.

 This plugin returns two scores:
 - pHaplo score gives the probability of a gene being haploinsufficient (deletion intolerant)
 - pTriplo score gives the probability of a gene being triploinsensitive (duplication intolerant)

 Pre-requisites:
 You need the compressed tsv file containing the dosage sensitivity score. The file
 Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz can be downloaded from here - 
 https://zenodo.org/record/6347673

 Options are passed to the plugin as key=value pairs:

 file   : (mandatory) compressed tsv file containing dosage sensitivity scores
 cover  : set value to 1 to report scores only if the variant covers the 
 affected gene completely (e.g. - a CNV that duplicates the gene); 0 by default
 
=cut

package DosageSensitivity;

use strict;
use warnings;

use Bio::EnsEMBL::VEP::Utils qw(get_compressed_filehandle);

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  my $param_hash = $self->params_to_hash();

  die "ERROR: file must be specified, such as file=/FULL_PATH_TO/dosage_sensitivity_scores.tsv.gz\n" unless defined $param_hash->{file};
  $self->{file} = $param_hash->{file};

  $self->{cover} = defined $param_hash->{cover} && $param_hash->{cover} eq '1' ? 1 : 0;

  $self->{header} = {
    "pHaplo" => "Probability of haploinsufficiency (deletion intolerance) of the affected gene",
    "pTriplo" => "Probability of triplosensitivity (duplication intolerance) of the affected gene"
  };

  return $self;
}

sub variant_feature_types {
    return ['VariationFeature', 'StructuralVariationFeature'];
}

sub feature_types {
  return ['Feature'];
}

sub get_header_info {
  my $self = shift;

  return $self->{header}
}

sub _process_DS_file {
  my ($self, $file) = @_;

  # get file handle
  my $fh = get_compressed_filehandle($file);

  my $dosage_sensitivity_matrix = {};
  my ($gene, $pHaplo, $pTriplo);
  while(<$fh>){
    next if $_[0] eq '#';

    ($gene, $pHaplo, $pTriplo) = split /\t/;
    $dosage_sensitivity_matrix->{$gene}->{pHaplo} = $pHaplo;
    $dosage_sensitivity_matrix->{$gene}->{pTriplo} = $pTriplo;
  }

  return $dosage_sensitivity_matrix;
}

sub _variant_cover_gene {
  my ($self, $tva) = @_;

  my $vf = undef;
  if (ref($tva) eq "Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele"){
    $vf = $tva->structural_variation_feature;
  }
  else {
    $vf = $tva->variation_feature;
  }

  my $gene = $tva->transcript->get_Gene;

  return 0 unless (defined $vf && defined $gene);

  return $vf->start <= $gene->start && $vf->end >= $gene->end;
}

sub run {
  my ($self, $tva) = @_;
  
  my $gene_name = $tva->transcript->get_Gene->external_name;
  return {} unless defined $gene_name;
  
  if ($self->{cover}) {
    return {} unless $self->_variant_cover_gene($tva);
  }

  my $dosage_sensitivity_matrix = $self->_process_DS_file($self->{file});
  return $dosage_sensitivity_matrix->{$gene_name};
}

1;
