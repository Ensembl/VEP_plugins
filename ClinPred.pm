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
 ClinPred

=head1 SYNOPSIS
 
 mv ClinPred.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin ClinPred,/path/to/ClinPred/data.txt.gz

=head1 DESCRIPTION
 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that adds pre-calculated scores from ClinPred.
 ClinPred is a prediction tool to identify disease-relevant nonsynonymous variants.

 This Plugin is only available for GRCh37
 
 Please cite the ClinPred publication alongside the VEP if you use this resource:
 https://www.sciencedirect.com/science/article/pii/S0002929718302714?via%3Dihub 

 ClinPred scores can be downloaded from 
 https://sites.google.com/site/clinpred/download
 
 The following steps are neccessary to tabix the ClinPred.txt.gz file before running the plugin:
 gzip -d ClinPred.txt.gz # to unzip the text file 
 cat ClinPred.txt | tr " " "\t" > ClinPred_tabbed.tsv # to change the file to a tabbed delimited file 
 sed '1s/.*/#&/'  ClinPred_tabbed.tsv > tabbed_ClinPred.tsv  # to add a # in the first line of the file 
 sed '1s/C/c/' tabbed_ClinPred.tsv > ClinPred_tabbed.tsv # to convert the Chr to chr 
 bgzip ClinPred_tabbed.tsv 

 tabix -f -s 1 -b 2 -e 2 ClinPred_tabbed.tsv.gz

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
  
  my $file = $self->params->[0];
  $self->add_file($file);
  
  my $assembly = $self->{config}->{assembly};

  if ($assembly ne "GRCh37") {
    die "Assembly is not GRCh37, ClinPred only works with GRCh37. \n";
  }
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
