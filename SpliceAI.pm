=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

 SpliceAI

=head1 SYNOPSIS

 mv SpliceAI.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin SpliceAI,/path/to/spliceai_scores.vcf.gz

=head1 DESCRIPTION

 A VEP plugin that retrieves pre-calculated annotations from SpliceAI.
 SpliceAI is a deep neural network, developed in partnership with Illumina, 
 that predicts splice junctions from an arbitrary pre-mRNA transcript sequence.
 It is available for both GRCh37 and GRCh38.

 More information can be found at:
 https://pypi.org/project/spliceai/

 Please cite the SpliceAI publication alongside VEP if you use this resource:
 https://www.ncbi.nlm.nih.gov/pubmed/30661751

 Output: 
 The output includes the alt allele, gene symbol, delta scores (DS) and delta positions (DP) 
 for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). 
 The output format is: ALLELE:SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL


 The following steps are necessary before running this plugin:

 The files with the annotations for all possible substitutions (snv), 1 base insertions 
 and 1-4 base deletions (indel) within genes are available here:
 https://basespace.illumina.com/analyses/194103939/files?projectId=66029966

 GRCh37:
 tabix -p vcf spliceai_scores.raw.snv.hg37.vcf.gz
 or
 tabix -p vcf spliceai_scores.raw.indel.hg37.vcf.gz

 GRCh38:
 tabix -p vcf spliceai_scores.raw.snv.hg38.vcf.gz
 or
 tabix -p vcf spliceai_scores.raw.indel.hg38.vcf.gz

 The plugin can then be run:
 ./vep -i variations.vcf --plugin SpliceAI,/path/to/spliceai_scores.raw.indel.hg38.vcf.gz


=cut

package SpliceAI;

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

  die("ERROR: SpliceAI scores input file not specified or found!\n") unless defined($self->params->[0]) && -e $self->params->[0];

  $self->{file} = $self->params->[0];

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {

  return{
     'SpliceAI_pred'  => 'SpliceAI predicted effect on splicing. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE:SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL',
  };
}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;
  my $chr = $vf->{chr};

  my $ref_allele;
  my $alt_allele;

  # convert to vcf format to compare the alleles
  if($vf->allele_string =~ /-/) {
    my $convert_to_vcf = $vf->to_VCF_record;
    $ref_allele = ${$convert_to_vcf}[3];
    $alt_allele = ${$convert_to_vcf}[4];
  }
  else { 
    my @alleles = split (/\//, $vf->allele_string, 2);
    $ref_allele = $alleles[0];
    $alt_allele = $alleles[1];
  }

  my $end = $vf->{end};
  my $start = $vf->{start};
  ($start, $end) = ($end, $start) if $start > $end;

  my @data = @{$self->get_data($chr, $start, $end)} if(defined $chr);

  return {} unless(@data);

  my $result_data = '';

  foreach my $data_value (@data) {

    if($data_value->{result}) {

      # Ref and alt alleles from SpliceAI file
      my $mm_ref = $data_value->{ref};
      my $mm_alt = $data_value->{alt};

      my @alt_allele_list;
      # Multiple alt alleles
      if ($alt_allele =~ /\//) {
        @alt_allele_list = split qr/\//, $alt_allele;
      }
      else {
        $alt_allele_list[0] = $alt_allele;
      }

      foreach my $alt (@alt_allele_list) {
        if($ref_allele eq $mm_ref && $alt eq $mm_alt) {
          $result_data = $result_data . "," . $data_value->{result};
        }
      }
    }
  }

  my %hash;
  my $result;

  if($result_data ne '') {
    $result_data =~ s/,//;
    $hash{'SpliceAI_pred'} = $result_data;
    $result = \%hash;
  }
  else {
    $result = {};
  }

  return $result;
}

# Parse data from SpliceAI file
sub parse_data {
  my ($self, $line) = @_;

  my ($chr, $start, $id, $ref, $alt, $qual, $filter, $info) = split /\t/, $line;

  $info =~ s/SpliceAI=//;
  my @info_splited = split (qr/\|/,$info, 2);
  my $allele = $info_splited[0];
  my $data = $info_splited[1];

  return {
    chr    => $chr,
    start  => $start,
    ref    => $ref,
    alt    => $alt,
    result => $allele . ":" . $data,
  };
}

1;