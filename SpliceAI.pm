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
 SpliceAI is a deep neural network, developed by Illumina, Inc 
 that predicts splice junctions from an arbitrary pre-mRNA transcript sequence.

 Delta score of a variant, defined as the maximum of (DS_AG, DS_AL, DS_DG, DS_DL), 
 ranges from 0 to 1 and can be interpreted as the probability of the variant being 
 splice-altering. The author-suggested cutoffs are:
   0.2 (high recall)
   0.5 (recommended)
   0.8 (high precision)

 This plugin is available for both GRCh37 and GRCh38.

 More information can be found at:
 https://pypi.org/project/spliceai/

 Please cite the SpliceAI publication alongside VEP if you use this resource:
 https://www.ncbi.nlm.nih.gov/pubmed/30661751

 Running options:
 (Option 1) By default, this plugin appends all scores from SpliceAI files.
 (Option 2) It can be specified a score cutoff between 0 and 1.

 Output: 
 The output includes the alt allele, gene symbol, delta scores (DS) and delta positions (DP) 
 for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). 
 The output format is: ALLELE:SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
 For VCF output files '|' is replaced by '&'.

 If plugin is run with option 2, the output also contains a flag: 'CUTOFF_PASS' if delta score
 passes the cutoff, 'CUTOFF_NO_PASS' otherwise. 

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
 ./vep -i variations.vcf --plugin SpliceAI,/path/to/spliceai_scores.raw.indel.hg38.vcf.gz,0.5


=cut

package SpliceAI;

use strict;
use warnings;
use List::Util qw(max);

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
  if(defined($self->params->[1])) {
    my $cutoff = $self->params->[1];
    if($cutoff < 0 || $cutoff > 1) {
      die("ERROR: Cutoff score must be between 0 and 1!\n");
    }
    $self->{cutoff} = $cutoff;
  }

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  my %header;
  $header{'SpliceAI_pred'} = 'SpliceAI predicted effect on splicing. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE:SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL';

  if($self->{cutoff}) {
    $header{'SpliceAI_cutoff'} = 'Flag if delta score pass the cutoff (CUTOFF_PASS) or if it does not (CUTOFF_NO_PASS)'; 
  }

  return \%header;
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
  my $result_flag;

  foreach my $data_value (@data) {

    if($data_value->{result}) {

      # Ref and alt alleles from SpliceAI file
      my $splice_ref = $data_value->{ref};
      my $splice_alt = $data_value->{alt};

      my @alt_allele_list;
      # Multiple alt alleles
      if ($alt_allele =~ /\//) {
        @alt_allele_list = split qr/\//, $alt_allele;
      }
      else {
        $alt_allele_list[0] = $alt_allele;
      }

      if($ref_allele eq $splice_ref) {
        foreach my $alt (@alt_allele_list) {
          if($alt eq $splice_alt) {
            $result_data = $result_data . "," . $data_value->{result};

            # Add a flag if cutoff is used
            if($self->{cutoff}) {
              if($data_value->{info} >= $self->{cutoff}) {
                $result_flag = 'CUTOFF_PASS';
              }
              else {
                $result_flag = 'CUTOFF_NO_PASS';
              }
            }
          }
        }
      }
    }
  }

  my %hash;
  my $result;

  if($result_data ne '') {
    $result_data =~ s/,//;
    $hash{'SpliceAI_pred'} = $result_data;
    if($self->{cutoff}) {
      $hash{'SpliceAI_cutoff'} = $result_flag;
    }
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

  my $max_score;
  if($self->{cutoff}){
    my @scores = split (qr/\|/,$data);

    my @scores_list;
    push @scores_list, $scores[1];
    push @scores_list, $scores[2];
    push @scores_list, $scores[3];
    push @scores_list, $scores[4];

    $max_score = max(@scores_list);
  }

  return {
    chr    => $chr,
    start  => $start,
    ref    => $ref,
    alt    => $alt,
    info   => $max_score,
    result => $allele . ":" . $data,
  };
}

1;