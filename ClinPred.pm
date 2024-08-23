=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
 ClinPred

=head1 SYNOPSIS
 
 mv ClinPred.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin ClinPred,file=ClinPred_tabbed.tsv.gz


=head1 DESCRIPTION
 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that adds pre-calculated scores from ClinPred.
 ClinPred is a prediction tool to identify disease-relevant nonsynonymous variants.
 
 Please cite the ClinPred publication alongside the VEP if you use this resource:
 https://www.sciencedirect.com/science/article/pii/S0002929718302714

 ClinPred scores can be downloaded from 
 https://sites.google.com/site/clinpred/download
 
 The following steps are neccessary to tabix the ClinPred.txt.gz file before running the plugin:

 For GRCh37:
 gzip -d ClinPred.txt.gz # to unzip the text file 
 cat ClinPred.txt | tr " " "\t" > ClinPred_tabbed.tsv # change to tab-delimited file 
 sed -i '1s/.*/#&/' ClinPred_tabbed.tsv  # comment the first line
 sed -i '1s/Chr/chr/' ClinPred_tabbed.tsv # convert Chr to chr
 bgzip ClinPred_tabbed.tsv
 tabix -f -s 1 -b 2 -e 2 ClinPred_tabbed.tsv.gz

 For GRCh38:
 gzip -d ClinPred_hg38.txt.gz # unzip the text file 
 awk '($2 == "Start" || $2 ~ /^[0-9]+$/){print $0}' ClinPred_hg38.txt > "ClinPred_hg38_tabbed.tsv" # remove problematic lines
 sed -i '1s/.*/#&/' ClinPred_hg38_tabbed.tsv # comment the first line
 sed -i '1s/Chr/chr/' ClinPred_hg38_tabbed.tsv # convert Chr to chr
 { head -n 1 ClinPred_hg38_tabbed.tsv; tail -n +2 ClinPred_hg38_tabbed.tsv | sort -k1,1V -k2,2V; } > ClinPred_hg38_sorted_tabbed.tsv # sort file by chromosome and position
 bgzip ClinPred_hg38_sorted_tabbed.tsv
 tabix -f -s 1 -b 2 -e 2 ClinPred_hg38_sorted_tabbed.tsv.gz

 The tabix utility must be installed in your path to use this plugin.
 Check https://github.com/samtools/htslib.git for instructions.


=cut

package ClinPred;
 
use strict;
use warnings; 

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my %INCLUDE_SO = map {$_ => 1} qw(missense_variant stop_lost stop_gained start_lost);
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

  $self->add_file($file);
  
  my $assembly = $self->{config}->{assembly};

  return $self;
}

sub feature_types {
  return ['Transcript'];

}

sub get_header_info {
  return { ClinPred => "Prediction tool to identify disease-relevant nonsynonymous single-nucleotide variants" };
}
 
sub run{
  my ($self, $tva) = @_;
  return {} unless grep {$INCLUDE_SO{$_->SO_term}} @{$tva->get_all_OverlapConsequences};
  
  my $vf = $tva->variation_feature;
  my $allele = $tva->variation_feature_seq;
  return {} unless $allele =~ /^[ACGT]$/;
  

  my ($vf_start, $vf_end) = ($vf->{start}, $vf->{end});
  ($vf_start, $vf_end) = ($vf_end, $vf_start) if ($vf_start > $vf_end);

  my ($res) = grep{
    $_->{alt} eq $allele &&
    $_->{start} == $vf_start &&
    $_->{end} == $vf_end
    } @{$self->get_data($vf->{chr}, $vf_start, $vf_end)};
    
  return $res ? $res->{result} : {};
}

sub parse_data {
  my ($self, $line) = @_;
  my ($chr, $start, $ref, $alt, $clinpred_score ) = split("\t", $line);
  return {
    start => $start,
    end => $start,
    alt => $alt,
    result => {
      ClinPred => $clinpred_score,
    }
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
