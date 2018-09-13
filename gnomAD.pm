=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

 Stephen Kazakoff <sh.kazakoff@gmail.com>
    
=cut

=head1 NAME

 gnomAD

=head1 SYNOPSIS

 mv gnomAD.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin gnomAD,/path/to/gnomad-public/r2.0.2/vcf/genomes,col1,col2
 ./vep -i variations.vcf --plugin gnomAD,/path/to/gnomad-public/r2.0.2/vcf/exomes,col1,col2

=head1 DESCRIPTION

 A VEP plugin that retrieves gnomAD annotation from either the genome
 or exome VCF files, available here:

 http://gnomad.broadinstitute.org/downloads

 The genome or exome VCF files must be downloaded and tabix indexed before
 using this plugin:

 > release="2.0.2"

 > genomes="https://storage.googleapis.com/gnomad-public/release/${release}/vcf/genomes"
 > wget -x "${genomes}"/gnomad.genomes.r${release}.sites.chr{{1..22},X}.vcf.bgz{,.tbi}

 > exomes="https://storage.googleapis.com/gnomad-public/release/${release}/vcf/exomes"
 > wget -x "${exomes}"/gnomad.exomes.r${release}.sites.vcf.bgz{,.tbi}

 When running the plugin you must list at least one column to retrieve from the
 gnomAD VCF files, specified as parameters to the plugin e.g.

 --plugin gnomAD,/path/to/gnomad-public/r2.0.2/vcf/genomes,AF_POPMAX

 You may include all columns with ALL; this fetches a large amount of data per
 variant!:

 --plugin gnomAD,/path/to/gnomad-public/r2.0.2/vcf/genomes,ALL

 The parent directory's basename is used to set the output field prefix. This
 is 'gnomADg' for genomes, 'gnomADe' for exomes, or else just 'gnomAD'.

 If you use this plugin, please see the terms and data information:

 http://gnomad.broadinstitute.org/terms

 The gnomAD VCF files are provided for GRCh37, but if you use GRCh38
 you may like to use the liftover files, available here:

 https://console.cloud.google.com/storage/browser/gnomad-public/release

 The tabix utility must be installed in your path to use this plugin.


=cut

package gnomAD;

use strict;
use warnings;

use File::Basename;

use Bio::EnsEMBL::Variation::Utils::VEP qw(parse_line);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  # test tabix
  `tabix --version` || die("ERROR: tabix does not seem to be in your path\n");

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();

  my $params = $self->params;

  my $dir = shift @$params;
  die("ERROR: gnomAD directory not specified\n") unless $dir;
  die("ERROR: gnomAD directory not found\n") unless -d $dir;

  # use the parent directory's basename to set a column prefix
  my $base = basename($dir);
  my $prefix = 'gnomAD' . ($base eq 'genomes' || $base eq 'exomes' ? substr($base, 0, 1) : '');

  # add any VCF files to our list of inputs
  opendir (my $fh, $dir) or die $!;
  for (readdir $fh) {
    $self->add_file("$dir/$_") if /\.vcf\.b?gz$/;
  }
  closedir $fh;

  my @files = @{ $self->files() };

  die("ERROR: Could not find any $prefix VCF files\n") unless @files;

  my %all_headers;

  for my $file (@files) {
    my $info_fields = $self->parse_info_fields($file);
    my @info_field_headers = keys %$info_fields;

    die("ERROR: Could not read $prefix column headers from $file\n") unless scalar @info_field_headers > 1;

    @all_headers{@info_field_headers} = @{$info_fields}{@info_field_headers};
  }

  my $available_columns = join(',', sort keys %all_headers);

  my %required_headers = map { $_ => undef } @$params;

  if (exists $required_headers{CSQ}) {
    die("ERROR: Will not retrieve redundant CSQ annotations. Available $prefix columns are:\n$available_columns\n");
  }
  else {
    delete $all_headers{CSQ};
  }

  my %return_headers;

  if (exists $required_headers{ALL}) {
    delete $required_headers{ALL};

    %return_headers = %all_headers;
  }

  for (keys %required_headers) {
    if (defined(my $desc = $all_headers{$_})) {
      $return_headers{$_} = $desc;
    }
    else {
      die("ERROR: $_ column not found. Available $prefix columns are:\n$available_columns\n");
    }
  }

  unless (keys %return_headers) {
    die("ERROR: No columns selected to fetch. Available $prefix columns are:\n$available_columns\n");
  }

  $self->{prefix} = $prefix;
  $self->{headers} = \%return_headers;

  return $self;
}

sub parse_info_fields {
  my ($self, $file) = @_;

  my %headers;

  open(my $fh, '-|', "tabix -H $file") or die $!;

  while (<$fh>) {
    chomp;

    next unless /
      \#\#INFO=
      <
        ID=(?<id>.*?),
        Number=(?<number>.*?),
        Type=(?<type>.*?),
        Description="(?<description>.*?)"
      >
    /x;

    $headers{$+{id}} = $+{description};
  }

  close $fh;

  return \%headers;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my $prefix = $self->{prefix};
  my %headers = %{ $self->{headers} };

  my %header_info;

  my @keys = keys %headers;
  my @prefixed_keys = map { join('_', $prefix, $_) } @keys;

  @header_info{@prefixed_keys} = @headers{@keys};

  return \%header_info;
}

sub run {
  my ($self, $vfoa) = @_;

  my $vf = $vfoa->variation_feature;

  my ($res) = grep {
    $_->{start} eq $vf->{start} &&
    $_->{end} eq $vf->{end} &&
    $_->{alt} eq $vfoa->variation_feature_seq
  }
  @{ $self->get_data($vf->{chr}, $vf->{start} - 1, $vf->{end}) };

  return $res ? $res->{result} : {};
}

sub parse_data {
  my ($self, $line) = @_;

  my $prefix = $self->{prefix};
  my ($vf) = @{ parse_line({format => 'vcf', minimal => 1}, $line) };
  my ($ref, $alt) = split(/\//, $vf->allele_string);

  my %result = map {
    my ($k, $v) = split(/=/, $_, 2);
    join('_', $prefix, $k) => $v
  }
  grep {
    m/=/
  }
  split(/;/, (split(/\t/, $line))[-1]);

  return {
    chr => $vf->{chr},
    start => $vf->{start},
    end => $vf->{end},
    ref => $ref,
    alt => $alt,
    result => \%result,
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;

