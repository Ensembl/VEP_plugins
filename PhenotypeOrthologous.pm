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
 ./vep -i variations.vcf --plugin PhenotypeOrthologous,file=rat_human.gff3.gz,model=rat

=head1 DESCRIPTION

 A VEP plugin that retrieves phenotype information associated with orthologous genes from model organisms.

 The plugin presently only supports two model organisms, rat and mouse. 

 The PhenotypeOrthologous file can be downloaded from https://ftp.ensembl.org/pub/current_variation/PhenotypeOrthologous

 The plugin can be run: 
   The file option is mandatory to run this plugin 

  ./vep -i variations.vcf --plugin PhenotypeOrthologous,file=rat_mouse_orthologous.gff3.gz
 
  To return only results for rat :
    ./vep -i variations.vcf --plugin PhenotypeOrthologous,file=rat_mouse_orthologous.gff3.gz,model=rat

    ./vep -i variations.vcf --plugin PhenotypeOrthologous,file=rat_mouse_orthologous.gff3.gz,model=rattus_norvegicus
  To return only results for mouse:
    ./vep -i variations.vcf --plugin PhenotypeOrthologous,file=rat_mouse_orthologous.gff3.gz,model=mouse

    ./vep -i variations.vcf --plugin PhenotypeOrthologous,file=rat_mouse_orthologous.gff3.gz,model=mus_musculus

 The tabix utility must be installed in your path to use this plugin.
 Check https://github.com/samtools/htslib.git for instructions.

=cut 

package PhenotypeOrthologous;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);
 
  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();
    
  my $params = $self->params_to_hash();

  my $file = $params->{file};  
  my $model = $params->{model};

  die "File needs to be specified to run the PhenotypeOrthologus plugin. \n" if  (!$file);
  $self->add_file($file);

  if ( defined($model) && $model ne "rattus_norvegicus" && $model ne "rat" && $model ne "mouse" && $model ne "mus_musculus") {
    die "PhenotypeOrthologous plugin only works for rat and mouse \n"
  }
  
  if (defined($model)) {
    $model = "rat" if $model eq "rattus_norvegicus";
    $model = "mouse" if $model eq "mus_musculus";
  }
  
  $self->{model} = $model;
  $self->{params} = $params;


  return $self;

}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;
  
  my %header;


  if (keys(%{$self->{params}}) == 1)  {
    $header{"PhenotypeOrthologous_RatOrthologous_geneid"} = "PhenotypeOrthologous RatGene associated with Rat";
    $header{"PhenotypeOrthologous_RatOrthologous_phenotype"} = "PhenotypeOrthologous RatPhenotypes associated with orthologous genes in Rat";
    $header{"PhenotypeOrthologous_MouseOrthologous_geneid"} = "PhenotypeOrthologous MouseGene associated with Mouse";
    $header{"PhenotypeOrthologous_MouseOrthologous_phenotype"} = "PhenotypeOrthologous MousePhenotypes associated with orthologous genes in Mouse";
  }

  if (defined($self->{model}) && $self->{model} eq "rat" ) {
    $header{"PhenotypeOrthologous_RatOrthologous_geneid"} = "PhenotypeOrthologous RatGene associated with Rat";
    $header{"PhenotypeOrthologous_RatOrthologous_phenotype"} = "PhenotypeOrthologous RatPhenotypes associated with orthologous genes in Rat";
  }

  if (defined($self->{model}) && $self->{model} eq "mouse" ) {
    $header{"PhenotypeOrthologous_MouseOrthologous_geneid"} = "PhenotypeOrthologous MouseGene associated with Mouse";
    $header{"PhenotypeOrthologous_MouseOrthologous_phenotype"} = "PhenotypeOrthologous MousePhenotypes associated with orthologous genes in Mouse";
  }

  return \%header;

}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;

  my $transcript = $tva->transcript;
  my $gene_id =  $transcript->{_gene}->stable_id;

  my ($res) = grep {
    $_->{gene_id} eq $gene_id
  } @{$self->get_data($vf->{chr}, $vf->{start}, $vf->{end})};
 
  return $res ? $res->{result} : {};

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
  
  if (!$self->{model}) {
    return {
      start => $s,
      end => $e,
      gene_id => $data_fields{"gene_id"},
      result => {
        PhenotypeOrthologous_RatOrthologous_geneid => $data_fields{"Rat_gene_id"},
        PhenotypeOrthologous_RatOrthologous_phenotype => $data_fields{"Rat_Orthologous_phenotype"},
        PhenotypeOrthologous_MouseOrthologous_geneid => $data_fields{"Mouse_gene_id"},
        PhenotypeOrthologous_MouseOrthologous_phenotype => $data_fields{"Mouse_Orthologous_phenotype"},
      }
    };
  }

  if ($self->{model} eq "rat") {
    return {
      start => $s,
      end => $e,
      gene_id => $data_fields{"gene_id"},
      result => {
        PhenotypeOrthologous_RatOrthologous_geneid => $data_fields{"Rat_gene_id"},
        PhenotypeOrthologous_RatOrthologous_phenotype => $data_fields{"Rat_Orthologous_phenotype"},
      }
    };
  }

  if ($self->{model} eq "mouse") {
    return {
      start => $s,
      end => $e,
      gene_id => $data_fields{"gene_id"},
      result => {
        PhenotypeOrthologous_MouseOrthologous_geneid => $data_fields{"Mouse_gene_id"},
        PhenotypeOrthologous_MouseOrthologous_phenotype => $data_fields{"Mouse_Orthologous_phenotype"},
      }
    };
  }


}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;