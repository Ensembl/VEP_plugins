=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2026] EMBL-European Bioinformatics Institute

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

Evo2VED - annotate variants with Evo2 Variant Effect Database scores

=head1 SYNOPSIS

 mv Evo2VED.pm ~/.vep/Plugins
 ./vep -i variations.vcf --assembly GRCh38 --plugin Evo2VED,file=/FULL_PATH/evo2ved.GRCh38.v1.0.0.scores.tsv.gz
 ./vep -i variations.vcf --assembly GRCh38 --plugin Evo2VED,file=/FULL_PATH/evo2ved.GRCh38.v1.0.0.scores.tsv.gz,Evo2_7B_AvgRC_Delta=1
 ./vep -i variations.vcf --assembly GRCh38 --plugin Evo2VED,file=/FULL_PATH/evo2ved.GRCh38.v1.0.0.scores.tsv.gz,all=1

=head1 DESCRIPTION

 An Ensembl VEP plugin that retrieves precomputed Evo2VED scores for substitutions,
 insertions and deletions from a tabix-indexed, bgzip-compressed TSV file. Scores
 are matched by position and REF/ALT alleles, including VEP's trimmed indel
 representation, and can annotate any feature type or intergenic variant.

 Evo2VED: https://dnallmved.org/
 Please cite Zhu R et al., Evo2VED: a database of DNA large language model-derived
 scores for variant interpretation. Nucleic Acids Research (2026), gkag1076:
 https://doi.org/10.1093/nar/gkag1076

 Download the GRCh38 TSV release and its matching index from
 https://dnallmved.org/downloads and keep both files in the same directory:

 curl -fL -O https://huggingface.co/datasets/huthvincent/Evo2VED-scores/resolve/main/v1.0.0/evo2ved.GRCh38.v1.0.0.scores.tsv.gz
 curl -fL -O https://huggingface.co/datasets/huthvincent/Evo2VED-scores/resolve/main/v1.0.0/evo2ved.GRCh38.v1.0.0.scores.tsv.gz.tbi

 No reformatting is needed. The tabix utility from HTSlib must be in your PATH
 (installation: https://www.htslib.org/download/). To recreate a missing index:

 tabix -s 1 -b 2 -e 2 evo2ved.GRCh38.v1.0.0.scores.tsv.gz

 The file must keep its header line, which names the columns:

 #CHROM POS REF ALT RSID Variant_ID Evo2_7B_NoRC_Delta Evo2_7B_AvgRC_Delta Evo2_7B_WeightedRC_Delta Evo2_40B_NoRC_Delta Evo2_40B_AvgRC_Delta Evo2_40B_WeightedRC_Delta

 The first four columns must be in that order. Score columns are found by name.
 Positions are 1-based. Variants absent from the file receive no annotation.

 The v1.0.0 release is for GRCh38 only and covers chromosomes 1-22. The plugin
 adds no annotation when VEP reports another assembly. With --gff or --gtf
 annotation and no cache, VEP cannot tell the assembly, so give
 --assembly GRCh38 there.

 Evo2_40B_AvgRC_Delta is reported by default. Enable any additional score using
 its exact column name followed by =1, or use all=1 to enable all six scores.
 An explicit column-name=0 disables that score, including the default or a score
 enabled by all=1. Boolean parameters also accept true/false and yes/no.
 The data file may alternatively be passed as a positional parameter.

 Each variant is looked up with its own tabix query. On dense, sorted input,
 expand_right=N (for example expand_right=100000) fetches and caches the next N
 bases with each query, and cache_size sets how many such regions are kept.

 Scores retain the original column names and signed values: Delta = ALT score
 minus REF score. Negative values indicate that ALT is less expected by the
 model; positive values indicate that ALT is more expected. Missing scores (.
 or NA) are omitted. The 7B and 40B models each provide forward-only (NoRC),
 reverse-complement-averaged (AvgRC), and weighted reverse-complement
 (WeightedRC) scores. These are model scores, not calibrated probabilities.

=cut

package Evo2VED;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my @SCORES = (
  ['Evo2_7B_NoRC_Delta',       '7B, forward strand only'],
  ['Evo2_7B_AvgRC_Delta',      '7B, reverse-complement averaged'],
  ['Evo2_7B_WeightedRC_Delta', '7B, weighted reverse complement'],
  ['Evo2_40B_NoRC_Delta',       '40B, forward strand only'],
  ['Evo2_40B_AvgRC_Delta',      '40B, reverse-complement averaged'],
  ['Evo2_40B_WeightedRC_Delta', '40B, weighted reverse complement'],
);
my $DEFAULT_SCORE = 'Evo2_40B_AvgRC_Delta';
my $ASSEMBLY = 'GRCh38';

sub new {
  my ($class, $config, @params) = @_;

  my %known = map { $_ => 1 } ('file', 'all', 'expand_left', 'expand_right', 'cache_size', map { $_->[0] } @SCORES);
  foreach my $param (grep { $known{$_} } @params) {
    die "ERROR: Evo2VED parameter '$param' needs a value, for example $param=1\n";
  }

  # Normalise a positional file before the base class parses key=value options.
  die "ERROR: Specify only one Evo2VED file\n"
    if (grep { index($_, '=') == -1 || /^file=/ } @params) > 1;
  @params = map { index($_, '=') == -1 ? "file=$_" : $_ } @params;
  my $self = $class->SUPER::new($config, @params);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();

  my $params = $self->params_to_hash();
  die "ERROR: No Evo2VED file specified\n" unless $params->{file};
  die "ERROR: Evo2VED scores require $ASSEMBLY\n" if _wrong_assembly($config->{assembly});

  foreach my $param (keys %$params) {
    die "ERROR: Unknown Evo2VED parameter '$param'\n" unless $known{$param};
  }

  my $all = exists($params->{all}) ? _enabled('all', $params->{all}) : 0;
  my @selected = grep {
    my $name = $_->[0];
    exists($params->{$name}) ? _enabled($name, $params->{$name}) : ($all || $name eq $DEFAULT_SCORE);
  } @SCORES;
  die "ERROR: No Evo2VED scores selected\n" unless @selected;

  $self->add_file($params->{file});

  my $names = _header_names($params->{file});
  my %column = map { $names->[$_] => $_ } 0 .. $#$names;
  $self->{width} = scalar @$names;
  $self->{scores} = [map {
    my ($name, $description) = @$_;
    die "ERROR: Evo2VED file has no $name column\n" unless exists $column{$name};
    +{name => $name, description => $description, column => $column{$name}};
  } @selected];

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;
  return {} if $self->_disabled;
  return {map {
    $_->{name} => "Evo2VED signed ALT minus REF score (Evo2 $_->{description})"
  } @{$self->{scores}}};
}

sub run {
  my ($self, $tva) = @_;
  return {} if $self->_disabled;

  my $vf = $tva->variation_feature;
  my ($start, $end) = ($vf->{start}, $vf->{end});
  my $ref = $vf->ref_allele_string;
  my $allele = $tva->variation_feature_seq;

  # With --minimal, VEP trims each ALT of a multi-allelic line separately and
  # restores the line's coordinates before plugins run. Read the untrimmed alleles.
  if ($vf->{minimised} && ($vf->{allele_string} // '') =~ m{/.+/}) {
    my @alleles = split m{/}, $vf->{allele_string};
    ($ref, $allele) = ($alleles[0], $alleles[$tva->allele_number // 0]);
  }

  return {} unless defined($start) && defined($end);
  return {} unless defined($ref) && $ref =~ /^(?:[ACGT]+|-)$/;
  return {} unless defined($allele) && $allele =~ /^(?:[ACGT]+|-)$/;

  # Include the VCF anchor preceding VEP's trimmed indels. Insertions have
  # start = end + 1, so this also gives them a valid tabix query interval.
  my $query_start = $start > 1 ? $start - 1 : 1;
  my $query_end = $end > $query_start ? $end : $query_start;
  foreach my $record (@{$self->get_data($vf->{chr}, $query_start, $query_end)}) {
    my $matches = get_matched_variant_alleles(
      {ref => $ref, alts => [$allele], pos => $start, strand => $vf->strand},
      {ref => $record->{ref}, alts => [$record->{alt}], pos => $record->{start}}
    );
    return $record->{result} if @{$matches};
  }

  return {};
}

sub parse_data {
  my ($self, $line) = @_;
  return if $line =~ /^#/;
  $line =~ s/[\r\n]+$//;
  my @values = split /\t/, $line, -1;
  return unless @values >= $self->{width};

  my %result;
  foreach my $score (@{$self->{scores}}) {
    my $value = $values[$score->{column}];
    next if $value eq '' || $value eq '.' || $value eq 'NA';
    $result{$score->{name}} = $value;
  }

  return {
    chr => $values[0], start => $values[1], end => $values[1],
    ref => $values[2], alt => $values[3], result => \%result,
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

sub _wrong_assembly {
  my $assembly = shift;
  return defined($assembly) && $assembly ne $ASSEMBLY;
}

# VEP resolves an omitted --assembly from the cache after plugins are constructed.
# It asks for headers before it forks, so forked runs warn once.
sub _disabled {
  my $self = shift;
  if (!$self->{disabled} && _wrong_assembly($self->{config}->{assembly})) {
    warn "WARNING: Evo2VED scores require $ASSEMBLY; no Evo2VED annotation will be added\n";
    $self->{disabled} = 1;
  }
  return $self->{disabled};
}

sub _header_names {
  my $file = shift;

  open(my $fh, '-|', 'tabix', '-H', $file) or die "ERROR: Could not run tabix: $!\n";
  my @lines = <$fh>;
  close($fh) or die "ERROR: tabix could not read the header of $file\n";

  my ($header) = grep { /^#CHROM\tPOS\tREF\tALT\t/ } @lines;
  die "ERROR: $file has no '#CHROM POS REF ALT' header line\n" unless defined $header;

  $header =~ s/[\r\n]+$//;
  return [split /\t/, $header];
}

sub _enabled {
  my ($name, $value) = @_;
  return 1 if defined($value) && $value =~ /^(?:1|true|yes)$/i;
  return 0 if defined($value) && $value =~ /^(?:0|false|no)$/i;
  die "ERROR: Evo2VED parameter '$name' must be 0/1, false/true or no/yes\n";
}

1;
