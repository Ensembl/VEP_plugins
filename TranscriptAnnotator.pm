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

 TranscriptAnnotator

=head1 SYNOPSIS

 mv TranscriptAnnotator.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin TranscriptAnnotator,file=/path/to/file.txt.gz

=head1 DESCRIPTION

 A VEP plugin that annotates variant-transcript pairs based on a given file:

   --plugin TranscriptAnnotator,file=${HOME}/file.tsv.gz

 Example of a valid tab-separated annotation file:

 ```
 #Chrom  Pos       Ref  Alt  Transcript       SIFT_score  SIFT_pred    Comment
 11      436154    A    G    NM_001347882.2   0.03        Deleterious  Bad
 11      1887471   C    T    ENST00000421485  0.86        Tolerated    Good
 ```

 Please bgzip and tabix the file with commands such as:

   bgzip file.txt
   tabix -b2 -e2 file.txt.gz

 Options are passed to the plugin as key=value pairs:

   file   : (mandatory) Tabix-indexed file to parse. Must contain variant
            location (chromosome, position, reference allele, alternative allele)
            and transcript ID as the first 5 columns. Accepted transcript IDs
            include those from Ensembl and RefSeq.
   cols   : Colon-delimited list with names of the columns to append. Column
            names are based on the last header line. By default, all columns
            (except the first 5) are appended.
   prefix : String to prefix the name of appended columns (default: basename of
            the filename without extensions). Set to 0 to avoid any prefix.
   trim   : Trim whitespaces from both ends of each column (default: 1).

 The tabix and bgzip utilities must be installed in your path to read the
 tabix-indexed annotation file: check https://github.com/samtools/htslib.git for
 installation instructions.

=cut

package TranscriptAnnotator;

use strict;
use warnings;

use File::Basename;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub _basename_without_ext {
  my $file = basename shift;
  $file =~ s/\.[^.]*?(.gz)?$//g;
  return $file;
}

sub _trim_whitespaces {
  my $text = shift;
  $text =~ s/^\s+|\s+$//g;
  return $text;
};

sub _get_colnames {
  my $self = shift;

  # Open file header
  open IN, "tabix -f -h " . $self->{_files}[0] . " 1:1-1 |"
    or die "ERROR: cannot open tabix file for " . $self->{_files}[0];

  # Get last line
  my $last;
  $last = $_ while <IN>;
  $last =~ s/(^#|\n$)//g;
  close IN;

  # Parse column names from header
  my @cols = split /\t/, $last;
  @cols = map { _trim_whitespaces $_ } @cols if $self->{trim};

  # Check column names for variant location and transcript ID
  my @var_cols = ("chr", "pos", "ref", "alt", "tr");
  for my $i (0 .. scalar $#var_cols) {
    unless (grep { /^$var_cols[$i]/i } @cols) {
      warn sprintf "WARNING: Name of column %s ('%s') should start with '%s'\n",
                   $i, $cols[$i], $var_cols[$i];
    }
  }

  # Return other columns
  @cols = splice @cols, 5;
  return @cols;
}

sub _get_selected_cols {
  my $self = shift;
  my @colnames = @_;

  my @cols;
  my $param_hash = $self->params_to_hash();
  if ( $param_hash->{cols} ) {
    # Check if user-selected columns are valid
    my @invalid_cols;
    for my $col ( split /:/, $param_hash->{cols} ) {
      if (grep /^$col$/, @colnames) {
        push @cols, $col;
      } else {
        push @invalid_cols, $col;
      }
    }

    # Warn about invalid user-selected columns
    my $msg_valid_cols = "Valid column names: " . join(", ", @colnames) . "\n";
    if (not @cols) {
      die "ERROR: no valid columns selected. " . $msg_valid_cols;
    } elsif (@invalid_cols) {
      my $file = $self->{_files}[0];
      warn "WARNING: columns " . join(", ", @invalid_cols) .
           " not found in $file. " . $msg_valid_cols;
    }
  } else {
    @cols = @colnames;
  }
  return @cols;
}

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);
  $self->get_user_params();

  my $param_hash = $self->params_to_hash();
  $self->{trim} = defined($param_hash->{trim}) ? $param_hash->{trim} : 1;

  # Check file
  my $file = $param_hash->{file};
  die "\nERROR: No file specified\nTry using 'file=path/to/file.txt.gz'\n"
     unless defined($file);
  $self->add_file($file);

  # Prepare column names
  my @colnames = $self->_get_colnames();
  my @cols     = $self->_get_selected_cols(@colnames);

  # Add prefix to column names
  unless (defined $param_hash->{prefix} && $param_hash->{prefix} eq 0) {
    $self->{prefix} = ($param_hash->{prefix} || _basename_without_ext($file)) . "_";
    @colnames = map($self->{prefix} . $_, @colnames);
    @cols     = map($self->{prefix} . $_, @cols);
  }
  $self->{colnames} = \@colnames;
  $self->{cols} = \@cols;

  return $self;
}

sub feature_types {
  return ['Transcript'];
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

sub run {
  my ($self, $tva) = @_;
  my $tr = $tva->transcript;

  # Get transcript ID for Ensembl and RefSeq
  my @refseq = split(/,/, $tr->{_refseq}) unless $tr->{_refseq} eq '-';
  my @tr_id  = ( $tr->{stable_id}, @refseq );
  my $vf     = $tva->variation_feature;

  # Get allele
  my $alt_alleles = $tva->base_variation_feature->alt_alleles;
  my $ref_allele  = $vf->ref_allele_string;

  my @data = @{ $self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end}) };
  return {} unless(@data);

  foreach my $var (@data) {
    my $is_same_transcript = grep { $var->{transcript_id} eq $_ } @tr_id;

    my $matches = get_matched_variant_alleles(
      {
        ref    => $ref_allele,
        alts   => $alt_alleles,
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
        ref  => $var->{ref},
        alts => [$var->{alt}],
        pos  => $var->{start},
      }
    );
    return $var->{result} if $is_same_transcript && (@$matches);
  }
  return {};
}

sub parse_data {
  my ($self, $line) = @_;

  my @data = split /\t/, $line;
  @data = map(_trim_whitespaces($_), @data) if $self->{trim};
  my ($seqname, $pos, $ref, $alt, $transcript_id, @vals) = @data;

  my %res;
  @res{ @{ $self->{colnames} } } = @vals;

  # Filter by user-selected columns
  %res = map { $_ => $res{$_} } @{ $self->{cols} };

  return {
    seqname => $seqname,
    start => $pos,
    end => $pos,
    ref => $ref,
    alt => $alt,
    transcript_id => $transcript_id,
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
