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

PhenotypeOrthologous

=head1 SYNOPSIS

 mv PhenotypeOrthologous.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin PhenotypeOrthologous

=head1 DESCRIPTION

 A VEP plugin that retrieves phenotype information associated with orthologous genes from from model organisms

=cut 

package PhenotypeOrthologous;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();
    
  my $params = $self->params_to_hash();
  my $file;
   
  if (!keys %$params) {
    $file = $self->params->[0];
    $params->{file} = $file;
  } else {
    $file = $params->{file};  
  }

  die "File needs to be specified to run the PhenotypeOrthologus plugin. \n" if  (!$file);
  $self->add_file($file);
  
  return $self;

}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return { PhenotypeOrthologous => 'Phenotypes associated with orthologous genes in model organism'};
}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;

  my $transcript = $tva->transcript;
  my $gene_id =  $transcript->{_gene}->stable_id;

  my @data = @{ $self->get_data($vf->{chr}, $vf->{start}, $vf->{end}) };
  return {} unless(@data);

  foreach (@data) {
    return $_->{result} if $_->{gene_id} eq $gene_id;
  }

  return {};
}


sub parse_data {
  my ($self, $line) = @_;

  my ($c, $grc, $feat, $s, $e, $n, $str, $n2, $note) = split /\t/, $line;
  
  my @fields = split /;/, $note;
  my %data_fields;
  
  foreach my $field (@fields){
    my ($key, $value) = split /=/, $field;
    $data_fields{$key} = $value;
  }

  return {
    chr => $c,
    start => $s,
    end => $e,
    gene_id => $data_fields{"gene_id"},
    result => $data_fields{"Rat_Orthologous_phenotype"},
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;