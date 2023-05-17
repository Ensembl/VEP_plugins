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

 MaveDB

=head1 SYNOPSIS

 mv CAPICE.pm ~/.vep/Plugins
 
 # print only scores from MaveDB data (default)
 ./vep -i variations.vcf --plugin MaveDB,file=/full/path/to/data.csv.gz
 
 # print all columns from MaveDB data
 ./vep -i variations.vcf --plugin MaveDB,file=/full/path/to/data.csv.gz,cols=all

=head1 DESCRIPTION

 A VEP plugin that retrieves data from MaveDB, a database that contains
 multiplex assays of variant effect, including deep mutational scans and
 massively parallel report assays.
 
 Options are passed to the plugin as key=value pairs:
   file             : (mandatory) Tabix-indexed MaveDB file
   cols             : (optional) Colon-separated columns to print from MaveDB
                      files; if set to 'all', all columns are printed
                      (default: 'urn:score')
   transcript_match : Only print data if (Ensembl or RefSeq) transcript
                      identifiers match to those from MaveDB (default: 0)
 
 Please cite the MaveDB publication alongside the VEP if you use this resource:
 https://doi.org/10.1186/s13059-019-1845-6
 
 The tabix utility must be installed in your path to use this plugin.

=cut

package MaveDB;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub _get_colnames {
  my $self = shift;

  # Open file header
  open IN, "tabix -f -h " . $self->{_files}[0] . " 1:1-1 |"
    or die "ERROR: cannot open tabix file for " . $self->{_files}[0];

  # Get last line from header
  my $last;
  $last = $_ while <IN>;
  $last =~ s/(^#|\n$)//g;
  close IN;

  # Parse column names from header
  my @cols = split /\t/, $last;
  @cols = splice @cols, 5; # first columns only identify the variant

  # Prefix all column names with "MaveDB_"
  @cols = map { $_ =~ /^MaveDB_/ ? $_ : "MaveDB_" . $_ } @cols;
  return \@cols;
}

sub _parse_colnames {
  my $self = shift;
  my $param_hash = $self->params_to_hash();

  # Parse file columns
  $self->{colnames} = $self->_get_colnames();
  my $cols = $param_hash->{cols} || "MaveDB_urn:MaveDB_score";
  
  if ($cols eq "all") {
    $self->{cols} = $self->{colnames};
  } else {
    my @cols = split(/:/, $cols);

    # Prefix all column names with "MaveDB_"
    @cols = map { $_ =~ /^MaveDB_/ ? $_ : "MaveDB_" . $_ } @cols;
    $self->{cols} = \@cols;

    # Check validity of all columns
    my @invalid_cols = grep { !($_ ~~ $self->{colnames}) } @cols;
    die "\n  ERROR: The following columns were not found in file header: ",
      join(", ", @invalid_cols), "\n" if @invalid_cols;
  }
}

sub new {
  my $class = shift;  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();

  my $param_hash = $self->params_to_hash();
  $self->{transcript_match} = $param_hash->{transcript_match} || 0;

  # Check file
  my $file = $param_hash->{file};
  die "\n  ERROR: No file specified\nTry using 'file=path/to/file.csv.gz'\n"
     unless defined($file);
  $self->add_file($file);

  $self->_parse_colnames();

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;
  my %header;
  my @keys = @{ $self->{colnames} };

  my @vals = map { "column from " . $self->{_files}[0] } @keys;
  @header{ @keys } = @vals;

  # Filter by user-selected columns
  %header = map { $_ => $header{$_} } @{ $self->{cols} };

  return \%header;
}

sub _aminoacid_changes_match {
  my ($self, $tva, $protein_var) = @_;

  my $aa_ref = substr($protein_var, 0, 1);
  my $aa_alt = substr($protein_var, -1);
  my $aa_pos = substr($protein_var, 1, -1);
  return 0 unless defined $aa_ref && defined $aa_alt && defined $aa_pos;

  my $vf_aa_ref = $tva->base_variation_feature_overlap->get_reference_TranscriptVariationAllele->peptide;
  my $vf_aa_alt = $tva->peptide;
  my $vf_aa_pos = $tva->base_variation_feature_overlap->translation_start;
  return 0 unless defined $vf_aa_ref && defined $vf_aa_alt && defined $vf_aa_pos;

  return $vf_aa_pos eq $aa_pos &&
         $vf_aa_ref eq $aa_ref &&
         $vf_aa_alt eq $aa_alt;
}

sub _transcripts_match {
  my ($self, $tva, $transcript) = @_;

  # Get transcript ID for Ensembl and RefSeq
  my $tr = $tva->transcript;
  my @refseq = split(/,/, $tr->{_refseq}) unless $tr->{_refseq} eq '-';
  my @tr_id  = ( $tr->{stable_id}, @refseq );

  return grep { $transcript eq $_ =~ s/\.[0-9]+//gr } @tr_id;
}

sub _join_results {
  my $all_results = shift;
  my $res = shift;

  # Create array of results per key
  for (keys %$res) {
    $all_results->{$_} = [] if !$all_results->{$_};
    push(@{ $all_results->{$_} }, $res->{$_} || "NA");
  }
  return $all_results;
}

sub run {
  my ($self, $tva) = @_;
  
  my $vf = $tva->variation_feature;

  # Get allele
  my $allele = $tva->base_variation_feature->alt_alleles;

  my @data = @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})};

  my $all_results = {};
  foreach (@data) {
    # Check if aminoacid changes match
    my $protein_var = $_->{result}->{MaveDB_protein_variant};
    if ($protein_var) {
      next unless $self->_aminoacid_changes_match($tva, $protein_var);
    }

    # Check if transcripts match
    if ($self->{transcript_match}) {
      my $transcript = $_->{result}->{MaveDB_refseq};
      next unless $self->_transcripts_match($tva, $transcript);
    }

    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => $allele,
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
       ref  => $_->{ref},
       alts => [$_->{alt}],
       pos  => $_->{start},
      }
    );

    # Filter user-selected columns
    my %res = %{ $_->{result} };
    %res = map { $_ => $res{$_} } @{ $self->{cols} };

    $all_results = _join_results($all_results, \%res) if (@$matches);
  }
  return $all_results;
}

sub parse_data {
  my ($self, $line) = @_;
  my ($chrom, $start, $end, $ref, $alt, @vals) = split /\t/, $line;

  # VCF-like adjustment of mismatched substitutions for comparison with VEP
  if(length($alt) != length($ref)) {
    my $first_ref = substr($ref, 0, 1);
    my $first_alt = substr($alt, 0, 1);
    if ($first_ref eq $first_alt) {
      $start++;
      $ref = substr($ref, 1);
      $alt = substr($alt, 1);
      $ref ||= '-';
      $alt ||= '-';
    }
  }
  
  my %res;
  @res{ @{ $self->{colnames} } } = @vals;

  return {
    ref => $ref,
    alt => $alt,
    start => $start,
    result => \%res
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;