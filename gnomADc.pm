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
 ./vep -i variations.vcf --plugin gnomADc,/path/to/gnomADc.gz

=head1 DESCRIPTION

 A VEP plugin that retrieves gnomAD annotation from either the genome
 or exome coverage files, available here:

 https://gnomad.broadinstitute.org/downloads

 Or via the Google Cloud console:

 https://console.cloud.google.com/storage/browser/gnomad-public/release

 The coverage summary files must be processed and Tabix indexed before
 use by this plugin. Please select from the instructions below:

 # GRCh38 and gnomAD genomes:
 > genomes="https://storage.googleapis.com/gnomad-public/release/3.0/coverage/genomes"
 > genome_coverage_tsv="gnomad.genomes.r3.0.coverage.summary.tsv.bgz"
 > wget "${genomes}/${genome_coverage_tsv}"
 > zcat "${genome_coverage_tsv}" | sed -e '1s/^locus/#chrom\tpos/; s/:/\t/' | bgzip > gnomADc.gz
 > tabix -s 1 -b 2 -e 2 gnomADc.gz

 # GRCh37 and gnomAD genomes:
 > genomes="https://storage.googleapis.com/gnomad-public/release/2.1/coverage/genomes"
 > genome_coverage_tsv="gnomad.genomes.coverage.summary.tsv.bgz"
 > wget "${genomes}/${genome_coverage_tsv}"
 > zcat "${genome_coverage_tsv}" | sed -e '1s/^/#/' | bgzip > gnomADg.gz
 > tabix -s 1 -b 2 -e 2 gnomADg.gz

 # GRCh37 and gnomAD exomes:
 > exomes="https://storage.googleapis.com/gnomad-public/release/2.1/coverage/exomes"
 > exome_coverage_tsv="gnomad.exomes.coverage.summary.tsv.bgz"
 > wget "${exomes}/${exome_coverage_tsv}"
 > zcat "${exome_coverage_tsv}" | sed -e '1s/^/#/' | bgzip > gnomADe.gz
 > tabix -s 1 -b 2 -e 2 gnomADe.gz

 By default, the output field prefix is 'gnomAD'. However if the input file's
 basename is 'gnomADg' (genomes) or 'gnomADe' (exomes), then these values are
 used instead. This makes it possible to call the plugin twice and include
 both genome and exome coverage values in a single run. For example:

 ./vep -i variations.vcf --plugin gnomADc,/path/to/gnomADg.gz --plugin gnomADc,/path/to/gnomADe.gz

 This plugin also tries to be backwards compatible with older versions of the
 coverage summary files, including releases 2.0.1 and 2.0.2. These releases
 make available one coverage file per chromosome and these can be used "as-is"
 without requiring any preprocessing. To annotate against multiple tabix-indexed
 chromosome files, instead specify the path to the parent directory. For example:

 ./vep -i variations.vcf --plugin gnomADc,/path/to/gnomad-public/release/2.0.2/coverage/genomes

 When a directory path is supplied, only files immediately under this directory
 that have a '.txt.gz' extension will attempt to be loaded. By default, the
 output field prefix is simply 'gnomAD'. However if the parent directory is
 either 'genomes' or 'exomes', then the output field prefix will be 'gnomADg'
 or 'gnomADe', respectively.

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

  my $path = shift @{ $self->params };
  die("ERROR: gnomADc input file not specified\n") unless $path;
  die("ERROR: gnomADc input file not found\n") unless -e $path;

  my $prefix = 'gnomAD';

  if (-f $path) {

    $self->add_file($path);

    # use the file's basename to set a column prefix
    my $basename = fileparse($path, qr/\.[^.]*/);
    my %gnomad_basenames = map { $_ => 1 } qw( gnomADg gnomADe );

    $prefix = $basename if exists $gnomad_basenames{$basename};
  }
  elsif (-d $path) {

    opendir (my $fh, $path) or die $!;
    for (readdir $fh) {
      $self->add_file(File::Spec->catfile($path, $_)) if /\.txt\.gz$/;
    }
    closedir $fh;

    # use the parent directory's basename to set a column prefix
    my $basename = basename($path);
    my %gnomad_dirnames = map { $_ => 1 } qw( genomes exomes );

    $prefix .= substr($basename, 0, 1) if exists $gnomad_dirnames{$basename};
  }

  my @files = @{ $self->files() };

  die("ERROR: Could not find any $prefix coverage files\n") unless @files;

  $self->{prefix} = $prefix;

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my $prefix = $self->{prefix};
  my %header_info;

  for (qw(mean median)) {
    $header_info{ join('_', $prefix, $_, 'cov') } = "$_ coverage";
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
  my ($chr, $pos, @cov) = split /\t/, $line;

  my @keys = map {
    join('_', $prefix, $_, 'cov')
  }
  qw(mean median 1x 5x 10x 15x 20x 25x 30x 50x 100x);

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

