=head1 LICENSE
Copyright [2026] EMBL-European Bioinformatics Institute

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
 promoterAI

=head1 SYNOPSIS
 mv promoterAI.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin promoterAI,/path/to/promoterAI.tsv.gz

=head1 DESCRIPTION
 An Ensembl VEP plugin that adds promoterAI scores to promoter variants, predicting their impact on gene expression.

 Options are passed to the plugin as key=value pairs:
   file : (mandatory) Tabix-indexed file from Illumina PromoterAI (see below)
   cols : (optional) Colon-separated list of columns to return from the plugin
          file (default: "tss_pos:promoterAI"); use 'all' to print all data
   match_to : (optional) Feature type to match promoterAI scores to.
          One of ["transcript", "gene"] (default: "transcript");

 To download the promoterAI scores file to use with VEP (GRCh38 based),
   please follow the instructions found in the README at https://github.com/Illumina/PromoterAI.
   You need a valid license agreement as described in the README to obtain and use the promoterAI scores.

 Please cite the promoterAI publication alongside Ensembl VEP if you use this resource:
 https://www.science.org/doi/10.1126/science.ads7373

 Necessary before using the plugin:
   Do the following steps to index the annotations file before using the plugin:
    zcat promoterAI.tsv.gz | sed '1s/.*/#&/' | bgzip > promoterAI.tsv.bgz
    tabix -s 1 -b 2 -e 2 promoterAI.tsv.bgz

 You must have the Bio::DB::HTS module or the tabix utility must be installed
 in your path to use this plugin.
=cut

package promoterAI;

use strict;
use warnings;

use File::Basename;
use File::Spec;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub _plugin_name {
  return 'PromoterAI';
}

sub _prefix_cols {
  my $cols = shift;
  my $prefix = _plugin_name() . '_';
  my @prefixed_cols = map { $_ =~ /^$prefix/ ? $_ : $prefix . $_ } @$cols;
  return \@prefixed_cols;
}

sub _get_colnames {
  my $self = shift;

  # Open file header
  open IN, "tabix -H " . $self->{_files}[0] . " |"
    or die "ERROR: cannot open tabix file for " . $self->{_files}[0];

  # Get last line from header
  my $last;
  $last = $_ while <IN>;
  $last =~ s/(^#|\n$)//g;
  close IN;

  # Parse column names from header
  my @cols = split /\t/, $last;
  @cols = splice @cols, 4; # first columns only identify the variant
  return \@cols;
}

sub _parse_colnames {
  my $self = shift;
  my $cols = shift;

  # Parse file columns
  $self->{colnames} = $self->_get_colnames();

  if ($cols eq "all") {
    $self->{cols} = $self->{colnames};
  } else {
    my @cols = split(/:/, $cols);
    $self->{cols} = \@cols;

    # Check validity of all columns
    my @invalid_cols;
    for my $col (@{ $self->{cols} }) {
      push(@invalid_cols, $col) unless grep(/^$col$/, @{ $self->{colnames} });
    }

    die "\n  ERROR: The following columns were not found in file header: ",
      join(", ", @invalid_cols), "\n" if @invalid_cols;
  }

  # Rename the "promoterAI" column to "score"
  $self->{colnames} = [ map { $_ eq "promoterAI" ? "score" : $_ } @{ $self->{colnames} } ];
  $self->{cols}     = [ map { $_ eq "promoterAI" ? "score" : $_ } @{ $self->{cols} } ];

  # Prefix all column names
  $self->{colnames} = _prefix_cols $self->{colnames};
  $self->{cols}     = _prefix_cols $self->{cols};
}

sub _join_results {
  my $self = shift;
  my $all_results = shift;
  my $res = shift;

  if ($self->{config}->{output_format} eq 'json' || $self->{config}->{rest}) {
    # Group results for JSON and REST
    my $name = _plugin_name();
    $all_results->{$name} = [] unless defined $all_results->{$name};
    push(@{ $all_results->{$name} }, $res);
  } else {
    # Create array of results per key
    for (keys %$res) {
      $all_results->{$_} = [] if !$all_results->{$_};
      push(@{ $all_results->{$_} }, $res->{$_} || "NA");
    }
  }
  return $all_results;
}

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();

  my $param_hash = $self->params_to_hash();

  my $file = $param_hash->{file};
  die "\n  ERROR: No promoterAI file specified\nTry using 'file=path/to/data.tsv.bgz'\n"
     unless defined($file);
  $self->add_file($file);

  # Parse column names
  my $cols = $param_hash->{cols} || "tss_pos:promoterAI";
  $self->_parse_colnames($cols);

  # Store the feature matching level
  my $match_to = $param_hash->{match_to} || "transcript";
  if($match_to ne 'transcript' && $match_to ne 'gene') {
    die "\n  ERROR: match_to must be 'transcript' or 'gene'\n";
  }
  $self->{match_to} = $match_to;

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  my %header;
  my @keys = @{ $self->{colnames} };

  my $description = "column from " . basename $self->{_files}[0];
  my @vals = map { $description } @keys;
  @header{ @keys } = @vals;

  # Filter by user-selected columns
  %header = map { $_ => $header{$_} } @{ $self->{cols} };

  return \%header;
}

sub _filter_selected_cols {
  my $res = shift;
  my $cols = shift;

  my %tmp = %$res;
  %tmp = map { $_ => $res->{$_} } @$cols;
  return \%tmp;
}

sub run {
  my ($self, $vfoa) = @_;

  my $vf = $vfoa->variation_feature;

  (my $vf_chr = $vf->{chr}) =~ s/^chr//;
  my ($vf_start, $vf_end) = ($vf->{start}, $vf->{end});

  $vf_end = $vf_start if $vf_start > $vf_end;

  my @data = @{ $self->get_data($vf_chr, $vf_start -2, $vf_end) };

  return {} unless @data;

  # Filter the data to match the VariantFeatureOverlapAllele's...
  my $joined_results = {};
  my $vfoa_variation_feature_seq = $vfoa->variation_feature_seq;
  my $vfoa_transcript_id = $vfoa->transcript->stable_id;
  my $vfoa_gene_id = $vfoa->transcript->{_gene_stable_id};

  for my $data_candidate (@data) {
    # * Alternative allele sequence
    my $candidate_alt_allele = $data_candidate->{'alt'};

    my $seq_filtered_matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => [ $vfoa_variation_feature_seq ],
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
        ref  => $data_candidate->{ref},
        alts => [ $candidate_alt_allele ],
        pos  => $data_candidate->{start},
        strand => $data_candidate->{strand}
      }
    );
    if(!@$seq_filtered_matches){
      next;
    }

    # * feature ID
    if( $self->{match_to} eq "transcript" ){
      my $candidate_transcript_id = $data_candidate->{'transcript_id'};
      $candidate_transcript_id =~ s/\.[0-9]+//; # remove transcript version

      if($vfoa_transcript_id ne $candidate_transcript_id){
        next;
      }
    }
    elsif( $self->{match_to} eq "gene" ){
      my $candidate_gene_id = $data_candidate->{'gene_id'};
      $candidate_gene_id =~ s/\.[0-9]+//; # remove gene version

      if($vfoa_gene_id ne $candidate_gene_id){
        next;
      }
    }

    my $promoterAI_result = $data_candidate->{'result'};
    # Filter user-selected columns
    $promoterAI_result = _filter_selected_cols($promoterAI_result, $self->{cols});

    $joined_results = $self->_join_results($joined_results, $promoterAI_result);
  }

  return $joined_results;
}

sub parse_data {
  my ($self, $line) = @_;

  my $header = $self->{file_column};
  my ($chr, $pos, $ref, $alt, $gene, $gene_id, $transcript_id, $strand, $tss_position, $score) = split /\t/, $line;

  # VCF-like adjustment of mismatched substitutions for comparison with VEP
  if(length($alt) != length($ref)) {
    my $first_ref = substr($ref, 0, 1);
    my $first_alt = substr($alt, 0, 1);
    if ($first_ref eq $first_alt) {
      $pos++;
      $ref = substr($ref, 1) || "-";
      $alt = substr($alt, 1) || "-";
    }
  }

  my %result;

  @result{ @{ $self->{colnames} } } = ($gene, $gene_id, $transcript_id, $strand, $tss_position, $score);

  return {
    ref => $ref,
    alt => $alt,
    start => $pos,
    strand => $strand,
    gene_id => $gene_id,
    transcript_id => $transcript_id,
    result => \%result
  }
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;

