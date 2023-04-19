=head1 LICENSE
Copyright [2018-2020] QIMR Berghofer Medical Research Institute

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
 Stephen Kazakoff <Stephen.Kazakoff@qimrberghofer.edu.au>
    
=cut

=head1 NAME
 gnomADc

=head1 SYNOPSIS
 mv gnomADc.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin gnomADc,/path/to/gnomad.tsv.gz 

=head1 DESCRIPTION
 A VEP plugin that retrieves gnomAD annotation from either the genome
 or exome coverage files, available here:
   https://gnomad.broadinstitute.org/downloads

 To download the gnomad coverage file in TSV format:
  for Assembly GRCh37: 
   gnomad genomes:
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1/coverage/genomes/gnomad.genomes.coverage.summary.tsv.bgz --no-check-certificate
   gnomad exomes:
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1/coverage/exomes/gnomad.exomes.coverage.summary.tsv.bgz --no-check-certificate
  
  for Assembly GRCh38: 
   gnomad genomes: 
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz --no-check-certificate 

 Necessary before using the plugin
  for Assembly GRCh37:
   The following steps are necessary to tabix the gnomad genomes coverage file :
    gunzip -c gnomad.genomes.coverage.summary.tsv.bgz | sed '1s/.*/#&/' > gnomad.genomes.tabbed.tsv
    bgzip gnomad.genomes.tabbed.tsv
    tabix -s 1 -b 2 -e 2 gnomad.genomes.tabbed.tsv.gz
   
   The following steps are neccessary to tabix the gnomad exomes coverage file :
    gunzip -c gnomad.exomes.coverage.summary.tsv.bgz | sed '1s/.*/#&/' > gnomad.exomes.tabbed.tsv
    bgzip gnomad.exomes.tabbed.tsv
    tabix -s 1 -b 2 -e 2 gnomad.exomes.tabbed.tsv.gz
 
 for Assembly GRCh38:
   The following steps are necessary to tabix the gnomad genomes coverage file :
    gunzip -c gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz | sed '1s/.*/#&/' > gnomad.genomesv3.tabbed.tsv
    sed "1s/locus/chr\tpos/; s/:/\t/g" gnomad.genomesv3.tabbed.tsv > gnomad.ch.genomesv3.tabbed.tsv
    bgzip gnomad.ch.genomesv3.tabbed.tsv
    tabix -s 1 -b 2 -e 2 gnomad.ch.genomesv3.tabbed.tsv
 
 This plugin also tries to be backwards compatible with older versions of the
 coverage summary files, including releases 2.0.1 and 2.0.2. These releases
 provide one coverage file per chromosome and these can be used "as-is"
 without requiring any preprocessing. 
 
 If you use this plugin, please see the terms and data information:
   https://gnomad.broadinstitute.org/terms
 
 You must have the Bio::DB::HTS module or the tabix utility must be installed
 in your path to use this plugin. 
=cut

package gnomADc;

use strict;
use warnings;

use File::Basename;
use File::Spec;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();

  my $file = $self->params->[0];
  

  my $prefix = 'gnomAD';
  my $headers;
  $self->add_file($file);
  
  open FH, "tabix -fh $file 1:1-1 2>&1 | ";
  while(<FH>){
      next unless /^\#/;
      chomp;
      $_ =~ s/^\#//;
      $self->{headers} = [split];
  }
    
  close(FH);
  $headers = scalar @{$self->{headers}};

  die("ERROR: Could not find any $prefix coverage files\n") unless $file;


  $self->{prefix} = $prefix;
  $self->{file_column} = $headers;
  
  
  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my $prefix = $self->{prefix};
  my $header = $self->{file_column} if defined($self->{file_column});
  my %header_info;
  

  if (defined $header && $header == 14) {
    for (qw(mean median_approx total_DP)) {
      $header_info{ join('_', $prefix, $_, 'cov') } = "$_ coverage";
    }
  }
  else {
    for (qw(mean median)) {
      $header_info{ join('_', $prefix, $_, 'cov') } = "$_ coverage";
    }
  }
  for (qw(1x 5x 10x 15x 20x 25x 30x 50x 100x)) {
    $header_info{ join('_', $prefix, $_, 'cov') } = "Fraction of samples at $_ coverage";
  }
  
  return \%header_info;
}

sub run {
  my ($self, $vfoa) = @_;

  my $vf = $vfoa->variation_feature;

  (my $vf_chr = $vf->{chr}) =~ s/^chr//;
  my ($vf_start, $vf_end) = ($vf->{start}, $vf->{end});

  $vf_end = $vf_start if $vf_start > $vf_end;

  my @results = @{ $self->get_data($vf_chr, $vf_start, $vf_end) };

  return {} unless @results;

  my %sums;

  # sum the values across each position
  for my $result (@results) {
    $sums{$_} += $result->{$_} for keys %$result;
  }

  # take the average of each of the values
  my %avgs = map { $_ => sprintf("%.4f", $sums{$_} / scalar @results) } keys %sums;

  return \%avgs;
}

sub parse_data {
  my ($self, $line) = @_;

  my $prefix = $self->{prefix};
  my $header = $self->{file_column};
  my ($chr, $pos, @cov) = split /\t/, $line;
  my @keys;
  if ($header == 14 ){
    @keys = map {
      join('_', $prefix, $_, 'cov')
    }
    qw(mean median_approx total_DP 1x 5x 10x 15x 20x 25x 30x 50x 100x);
  }
  else{
    @keys = map {
      join('_', $prefix, $_, 'cov')
    }
    qw(mean median 1x 5x 10x 15x 20x 25x 30x 50x 100x);
  }
  
  my %result;

  @result{@keys} = @cov;

  return \%result;
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;

