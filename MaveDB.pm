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

 mv MaveDB.pm ~/.vep/Plugins
 
 # print only scores for single aminoacid changes from MaveDB data (default)
 ./vep -i variations.vcf --plugin MaveDB,file=/full/path/to/data.csv.gz

 # print all scores associated with the genetic variant
 ./vep -i variations.vcf --plugin MaveDB,file=/full/path/to/data.csv.gz,single_aminoacid_changes=0

 # print all columns from MaveDB data
 ./vep -i variations.vcf --plugin MaveDB,file=/full/path/to/data.csv.gz,cols=all

=head1 DESCRIPTION

 A VEP plugin that retrieves data from MaveDB (https://www.mavedb.org), a
 database that contains multiplex assays of variant effect, including deep
 mutational scans and massively parallel report assays.

 To run the MaveDB plugin, please download the following files containing
 MaveDB data for GRCh38 (we do not currently host data for other assemblies):
 - https://ftp.ensembl.org/pub/current_variation/MaveDB/MaveDB_variants.tsv.gz
 - https://ftp.ensembl.org/pub/current_variation/MaveDB/MaveDB_variants.tsv.gz.tbi

 Options are passed to the plugin as key=value pairs:
   file                     : (mandatory) Tabix-indexed MaveDB file
   cols                     : Colon-separated columns to print from MaveDB
                              files; if set to 'all', all columns are
                              printed (default: 'urn:score:nt:pro')
   single_aminoacid_changes : Return matches for single aminoacid changes only;
                              if disabled, return all matches associated with a
                              genetic variant (default: 1)
   transcript_match         : Return results only if (Ensembl or RefSeq)
                              transcript identifiers match (default: 1)
 
 Please cite the MaveDB publication alongside the VEP if you use this resource:
 https://doi.org/10.1186/s13059-019-1845-6
 
 The tabix utility must be installed in your path to use this plugin.

=cut

package MaveDB;

use strict;
use warnings;

use File::Basename;
use Bio::SeqUtils;
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
  my $cols = shift;

  # Parse file columns
  $self->{colnames} = $self->_get_colnames();
  
  if ($cols eq "all") {
    $self->{cols} = $self->{colnames};
  } else {
    my @cols = split(/:/, $cols);

    # Prefix all column names with "MaveDB_"
    @cols = map { $_ =~ /^MaveDB_/ ? $_ : "MaveDB_" . $_ } @cols;
    $self->{cols} = \@cols;

    # Check validity of all columns
    my @invalid_cols;
    for my $col (@{ $self->{cols} }) {
      push(@invalid_cols, $col) unless grep(/^$col$/, @{ $self->{colnames} });
    }

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
  my $tr_match = $param_hash->{transcript_match};
  $self->{transcript_match} = defined $tr_match ? $tr_match : 1;

  my $aa_changes = $param_hash->{single_aminoacid_changes};
  $self->{single_aa_changes} = defined $aa_changes ? $aa_changes : 1;

  # Check file
  my $file = $param_hash->{file};
  die "\n  ERROR: No file specified\nTry using 'file=path/to/file.csv.gz'\n"
     unless defined($file);
  $self->add_file($file);

  # Parse column names
  my $cols = $param_hash->{cols} || "urn:score:nt:pro";
  $self->_parse_colnames($cols);

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;
  my %header;
  my @keys = @{ $self->{colnames} };

  my @vals = map { "column from " . basename $self->{_files}[0] } @keys;
  @header{ @keys } = @vals;

  # Custom headers
  $header{"MaveDB_score"}    = "MaveDB score - see MaveDB for interpretation of scores; " . $header{"MaveDB_urn"};
  $header{"MaveDB_nt"}  = "MaveDB HGVS (nucleotide); " . $header{"MaveDB_urn"};
  $header{"MaveDB_pro"} = "MaveDB HGVS (protein); " . $header{"MaveDB_urn"};
  $header{"MaveDB_urn"} = "MaveDB database identifier; " . $header{"MaveDB_urn"};

  # Filter by user-selected columns
  %header = map { $_ => $header{$_} } @{ $self->{cols} };

  return \%header;
}

sub _aminoacid_changes_match {
  my ($self, $tva, $protein_var) = @_;

  my ($aa_ref, $aa_pos, $aa_alt) = $protein_var =~ /([A-Za-z]+)([0-9]+)([A-Za-z=*]+)/;
  return 0 unless defined $aa_ref && defined $aa_alt && defined $aa_pos;

  $aa_alt = $aa_ref if $aa_alt eq "="; # silent mutation
  $aa_alt = "Ter"   if $aa_alt eq "*"; # nonsense mutation

  # return no match if cannot fetch reference allele
  my $bvf = $tva->base_variation_feature_overlap;
  return 0 unless $bvf->can('get_reference_TranscriptVariationAllele');

  my $vf_aa_ref = $bvf->get_reference_TranscriptVariationAllele->peptide;
  my $vf_aa_alt = $tva->peptide;
  my $vf_aa_pos = $tva->base_variation_feature_overlap->translation_start;
  return 0 unless defined $vf_aa_ref && defined $vf_aa_alt && defined $vf_aa_pos;

  my $vf_aa3_ref = Bio::SeqUtils->seq3(
    Bio::PrimarySeq->new('-seq'=> $vf_aa_ref, '-alphabet'=>'protein'));
  my $vf_aa3_alt = Bio::SeqUtils->seq3(
    Bio::PrimarySeq->new('-seq'=> $vf_aa_alt, '-alphabet'=>'protein'));

  return ($vf_aa_pos eq $aa_pos) &&
         ($vf_aa_ref eq $aa_ref || $vf_aa3_ref eq $aa_ref) &&
         ($vf_aa_alt eq $aa_alt || $vf_aa3_alt eq $aa_alt);
}

sub _transcripts_match {
  my ($self, $tva, $transcript) = @_;

  # return no match if cannot fetch transcript from TVA
  return 0 unless $tva->can('transcript');

  # Get transcript ID for Ensembl and RefSeq
  my $tr = $tva->transcript;
  my @refseq = split(/,/, $tr->{_refseq}) if defined $tr->{_refseq} and $tr->{_refseq} ne '-';
  my @tr_id  = ( $tr->{stable_id}, @refseq );

  return grep { $transcript eq $_ =~ s/\.[0-9]+//gr } @tr_id;
}

sub _join_results {
  my $self = shift;
  my $all_results = shift;
  my $res = shift;

  # Filter user-selected columns
  my %tmp = %$res;
  %tmp = map { $_ => $res->{$_} } @{ $self->{cols} };
  $res = \%tmp;

  if ($self->{config}->{output_format} eq 'json' || $self->{config}->{rest}) {
    # Group results for JSON and REST
    $all_results->{"MaveDB"} = [] unless defined $all_results->{"MaveDB"};
    push(@{ $all_results->{"MaveDB"} }, $res);
  } else {
    # Create array of results per key
    for (keys %$res) {
      $all_results->{$_} = [] if !$all_results->{$_};
      push(@{ $all_results->{$_} }, $res->{$_} || "NA");
    }
  }
  return $all_results;
}

sub _is_single_aa_change {
  my $hgvsp = shift;
  return $hgvsp !~ /\;/;
}

sub run {
  my ($self, $tva) = @_;
  
  my $vf = $tva->variation_feature;

  # Get allele
  my $allele = $tva->base_variation_feature->alt_alleles;

  # Increment 1 to account for insertions
  my @data = @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end} + 1)};

  my $all_results = {};
  foreach (@data) {
    # Check if scores are associated with single aminoacid changes
    my $hgvsp = $_->{result}->{MaveDB_pro};
    if ($self->{single_aa_changes}) {
      next if defined $hgvsp && !_is_single_aa_change($hgvsp);
    }

    # Check if aminoacid changes match (if single aminoacid changes only)
    my $protein_var = $_->{result}->{MaveDB_protein_variant};
    $protein_var ||= $hgvsp if defined $hgvsp && _is_single_aa_change($hgvsp);
    if ($protein_var) {
      next unless $self->_aminoacid_changes_match($tva, $protein_var);
    }

    # Check if transcripts match
    if ($self->{transcript_match}) {
      my $transcript = $_->{result}->{MaveDB_refseq};
      next unless $transcript eq '' or $self->_transcripts_match($tva, $transcript);
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
    $all_results = $self->_join_results($all_results, $_->{result}) if @$matches;
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
