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

 AVADA

=head1 SYNOPSIS
./vep -i variations.vcf --plugin AVADA,file=path/to/file,feature_match_by=<gene_symbol/ensembl_id/refseq_id>

=head1 DESCRIPTION
wget http://bejerano.stanford.edu/AVADA/avada_v1.00_2016.vcf.gz
gzip -d avada_v1.00_2016.vcf.gz
bgzip -c avada_v1.00_2016.vcf > avada_v1.00_2016.vcf.gz
tabix -p vcf avada_v1.00_2016.vcf.gz
Output: 
 The output includes two columns  AVADA_FEATURE_ID; AVADA_PMID
=cut

package AVADA;

use strict;
use warnings;
use Data::Dumper;
use List::MoreUtils qw(uniq);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);



my $feature_match_by;

sub get_header_info {
  my $self = shift;
  my %header;

  $header{"AVADA_PMID"} = "PubMed ID evidence for the variant as sreported by AVADA";
  $header{"AVADA_FEATURE_ID"} = "Feature ID associated with variant as reported by AVADA" ;  
  return \%header
}

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();
  my $param_hash = $self->params_to_hash();
  $self->add_file($param_hash->{file}); 
  $feature_match_by = $param_hash->{feature_match_by} || "gene_symbol"; 
  if ($feature_match_by eq "refseq_id") {
    $self->{config}->{refseq} = 1;
    # DB mode will not work for RefSeq ID
    die("ERROR: Matching by RefSeq ID is not possible with --database option") if ($self->{config}->{database});
  }
  return $self;
}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;
  my $tv = $tva->transcript_variation;
  my $end = $vf->{end};
  my $start = $vf->{start};
  my $transcript = $tva->transcript;
  my %output;
  my @data;

  if ( $feature_match_by eq "gene_symbol"){
    @data = grep {
    $_->{AVADA_GENE_SYMBOL} eq $transcript->{_gene_symbol}
    }@{$self->get_data($vf->{chr}, $start, $end)};
    $output{"AVADA_FEATURE_ID"} = $transcript->{_gene_symbol};
  }
  elsif ( $feature_match_by eq "refseq_id" ){
    my $refseq_transcript = $transcript->{stable_id};
    my $refseq_protein = $transcript->translation->{stable_id};
    @data = grep {
      $_->{AVADA_REFSEQ_ID} eq $refseq_transcript
      }@{$self->get_data($vf->{chr}, $start, $end)};
    if (scalar @data) {
      $output{"AVADA_FEATURE_ID"} = $refseq_transcript;
    }
    else{
    @data = grep {
      $_->{AVADA_REFSEQ_ID} eq $refseq_protein
      }@{$self->get_data($vf->{chr}, $start, $end)};
      $output{"AVADA_FEATURE_ID"} = $refseq_protein;
    } 
  }
  elsif ( $feature_match_by eq "ensembl_id" )
  {
  @data = grep {
    $_->{AVADA_ENSEMBL_ID} eq $transcript->{_gene_stable_id}
    }@{$self->get_data($vf->{chr}, $start, $end)};
  $output{"AVADA_FEATURE_ID"} = $transcript->{_gene_stable_id};
  }
  return {} unless scalar @data;
  my $pmid_string ;
  my %seen;
  foreach my $data_value (uniq @data) {
    if ( exists $seen{$data_value->{AVADA_PMID}} )
    {
      if ( $seen{$data_value->{AVADA_PMID}} <= 1 )
      {
      $pmid_string = $pmid_string ? $pmid_string.",".$data_value->{AVADA_PMID} : $data_value->{AVADA_PMID} ;
      }
      else
      {
      $seen{$data_value->{AVADA_PMID}}++ ;
      }
    }
    else
    {
    $seen{$data_value->{AVADA_PMID}} = 1;
    }
  }
  
  $output{"AVADA_PMID"} = $pmid_string;
  return \%output;
}

sub parse_data {
  my ($self, $line) = @_;
  my ($chr, $start, $id, $ref, $alt, $x, $xx, $data) = split /\t/, $line;
  my ($pmid, $gene_ensembl_id, $gene_symbol, $gene_refseq_id);
  my @data_split = split /;/, $data;
  foreach my $value (@data_split){
      $pmid = $value if $value =~ /PMID/;
      $gene_ensembl_id = $value if $value =~ /ENSEMBL_ID/; 
      $gene_symbol = $value if $value =~ /GENE_SYMBOL/; 
      $gene_refseq_id = $value if $value =~ /REFSEQ_ID/; 

  }
  $pmid =~ s/PMID=//;
  $gene_ensembl_id =~ s/ENSEMBL_ID=//;
  $gene_symbol =~ s/GENE_SYMBOL=//;
  $gene_refseq_id =~ s/REFSEQ_ID=//;
  return {
    AVADA_PMID => $pmid,
    AVADA_ENSEMBL_ID => $gene_ensembl_id,
    AVADA_GENE_SYMBOL => $gene_symbol,
    AVADA_REFSEQ_ID => $gene_refseq_id
  };
}

1;
