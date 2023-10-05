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
 ./vep -i variations.vcf --plugin DosageSensitivity,file=/FULL_PATH_TO/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz,sv_only=1

=head1 DESCRIPTION


 
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

  die "ERROR: file is not specified which is a mandatory parameter\n" unless defined $param_hash->{file};
  $self->{file} = $param_hash->{file};

  $self->{header} = {
    "pHaplo" => "haploinsufficiency likelihood score",
    "pTriplo" => "triplosensitivity likelihood score"
  };

  return $self;
}

sub variant_feature_types {
    return ['VariationFeature', 'StructuralVariationFeature'];
}

sub feature_types {
  return ['Transcript'];
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

sub run {
  my ($self, $tva) = @_;
  
  my $gene_name = $tva->transcript->get_Gene->external_name;
  print($gene_name,"\n");
  return {} unless defined $gene_name;

  my $dosage_sensitivity_matrix = $self->_process_DS_file($self->{file});
  use Data::Dumper;
  print Dumper($dosage_sensitivity_matrix->{$gene_name}), "\n";
  return $dosage_sensitivity_matrix->{$gene_name};
}

1;
