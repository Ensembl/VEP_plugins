=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

 LOEUF

=head1 SYNOPSIS

 mv LOEUF.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin LOEUF,file=/path/to/loeuf/data.tsv.gz,match_by=gene
 ./vep -i variations.vcf --plugin LOEUF,file=/path/to/loeuf/data.tsv.gz,match_by=transcript

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the LOEUF scores to VEP output. LOEUF stands for the "loss-of-function 
 observed/expected upper bound fraction." 

 NB: The plugin currently does not add the score for downstream_gene_variant and upstream_gene_variant 

 Please cite the LOEUF publication alongside the VEP if you use this resource:
 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7334197/


 LOEUF scores can be downloaded from
 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7334197/bin/41586_2020_2308_MOESM4_ESM.zip

 For GRCh37:

 These files can be tabix-processed by:
 Unzip 41586_2020_2308_MOESM4_ESM.zip
 cd supplement
 zcat supplementary_dataset_11_full_constraint_metrics.tsv.gz | (head -n 1 && tail -n +2  | sort -t$'\t' -k 76,76 -k 77,77n ) | bgzip -c > loeuf_dataset.tsv.gz
 tabix -f -S 1 -s 76 -b 77 -e 78 loeuf_dataset.tsv.gz

 The tabix utility must be installed in your path to use this plugin.

=cut

package LOEUF;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

use Scalar::Util qw(looks_like_number);


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  my $param_hash = $self->params_to_hash();

  # Get file
  die("ERROR: LOEUF file not provided or not found!\n") unless defined($param_hash->{file}) && -e $param_hash->{file};
  $self->add_file($param_hash->{file});

  # Check match_by argument
  if(defined($param_hash->{match_by})) {
    my $match_by = $param_hash->{match_by};
    $self->{match_by} = $match_by;
  }

  else{
    die("ERROR: Argument 'match_by' is undefined");
  }

  # Check assembly
  my $assembly = $self->{config}->{assembly};
  if ($assembly ne "GRCh37") {
    die "Assembly is not GRCh37, LOEUF only works with GRCh37. \n";
  }

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return { LOEUF => 'Loss-of-function observed/expected upper bound fraction'};
}

sub run {
  my ($self, $tva) = @_;

  return {} if grep {$_->SO_term eq 'downstream_gene_variant' || $_->SO_term eq 'upstream_gene_variant'} @{$tva->get_all_OverlapConsequences};
  
  my $vf = $tva->variation_feature;
  my $transcript = $tva->transcript;
  my $end = $vf->{end};
  my $start = $vf->{start};
  ($start, $end) = ($end, $start) if $start > $end;
  
  if ($self->{match_by} eq 'transcript'){
    my ($res) = grep {
    $_->{transcript_id}  eq $transcript->stable_id;
    } @{$self->get_data($vf->{chr}, $start, $end)};
    return $res ? $res->{result} : {};
  }

  elsif ($self->{match_by} eq 'gene'){
    my @data = grep {
    $_->{gene_id} eq $transcript->{_gene_stable_id}
    } @{$self->get_data($vf->{chr}, $start, $end)};
    # Get max loeuf value for the gene
    my $max_loeuf_score = 0;
    foreach (@data){
      if (looks_like_number($_->{result}->{LOEUF})){
        if ($_->{result}->{LOEUF} > $max_loeuf_score){
          $max_loeuf_score = $_->{result}->{LOEUF};
        }
      }
    }
    return $max_loeuf_score ? {LOEUF => $max_loeuf_score} : {};
  }

  else{
    return {};
  }

}

sub parse_data {
  my ($self, $line) = @_;
  my @values = split /\t/, $line;
  my ($transcript_id, $oe_lof_upper, $gene_id, $chromosome, $start_position, $end_position ) = @values[1,30,64,75,76,77];
  return {
    gene_id => $gene_id,
    transcript_id => $transcript_id,
    chromosome => $chromosome,
    start => $start_position,
    end => $end_position,
    result => {
      LOEUF   => $oe_lof_upper,
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
