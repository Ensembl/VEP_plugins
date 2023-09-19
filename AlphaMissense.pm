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

 AlphaMissense

=head1 SYNOPSIS

 mv AlphaMissense.pm ~/.vep/Plugins

 # print AlphaMissense scores and predictions (default)
 ./vep -i variations.vcf --plugin AlphaMissense,file=/full/path/to/file.tsv.gz

 # print all AlphaMissense information
 ./vep -i variations.vcf --plugin AlphaMissense,file=/full/path/to/file.tsv.gz,cols=all

 # only report results for the transcripts in the AlphaMissense prediction
 ./vep -i variations.vcf --plugin AlphaMissense,file=/full/path/to/file.tsv.gz,transcript_match=1

=head1 DESCRIPTION

 This plugin for the Ensembl Variant Effect Predictor (VEP) annotates missense variants with the 
 pre-computed AlphaMissense pathogenicity scores. AlphaMissense is a deep learning model developed 
 by Google DeepMind that predicts the pathogenicity of single nucleotide missense variants. 

 This plugin will add two annotations per missense variant:

 'am_pathogenicity', a continuous score between 0 and 1 which can be interpreted as the predicted 
 probability of the variant being pathogenic. 

 'am_class' is the classification of the variant into one of three discrete categories:  
 'Likely pathogenic', 'Likely benign', or 'ambiguous'. These are derived using the following 
  thresholds of am_pathogenicity: 
    'Likely benign' if am_pathogenicity < 0.34; 
    'Likely pathogenic' if am_pathogenicity > 0.564; 
    'ambiguous' otherwise. 

 These thresholds were chosen to achieve 90% precision for both pathogenic and benign ClinVar variants. 
 Note that AlphaMissense was not trained on ClinVar variants. Variants labeled as 'ambiguous' should be 
 treated as 'unknown' or 'uncertain' according to AlphaMissense. 

 This plugin is available for both GRCh37 (hg19) and GRCh38 (hg38) genome builds.

 The prediction scores of AlphaMissense can be downloaded from 
 https://console.cloud.google.com/storage/browser/dm_alphamissense 
 (AlphaMissense Database Copyright (2023) DeepMind Technologies Limited).  Data contained within the 
 AlphaMissense Database is provided for non-commercial research use only under CC BY-NC-SA 4.0 license. 
 Use of the AlphaMissense Database is subject to Google Cloud Platform Terms of Service  

 Please cite the AlphaMissense publication alongside the VEP
 if you use this resource: https://doi.org/10.1126/science.adg7492

 Disclaimer: The AlphaMissense Database and other information provided on or linked to this site is 
 for theoretical modelling only, caution should be exercised in use. It is provided "as-is" without 
 any warranty of any kind, whether express or implied. For clarity, no warranty is given that use of 
 the information shall not infringe the rights of any third party (and this disclaimer takes precedence 
 over any contrary provisions in the Google Cloud Platform Terms of Service). The information provided 
 is not intended to be a substitute for professional medical advice, diagnosis, or treatment, and does 
 not constitute medical or other professional advice.

 Before running the plugin for the first time, you need to create a tabix index (requires tabix to be 
 installed).

 > tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.tsv.gz

 > tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg19.tsv.gz



 Options are passed to the plugin as key=value pairs:
   file             : (mandatory) Tabix-indexed AlphaMissense data
   cols             : (optional) Colon-separated columns to print from
                      AlphaMissense data; if set to 'all', all columns are printed
                      (default: Missense_pathogenicity:Missense_class)
   transcript_match : Only print data if transcript identifiers match those from
                      AlphaMissense data (default: 0)

 AlphaMissense predictions are matched to input data by genomic location and protein change by default.


=cut

package AlphaMissense;

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
  @cols = splice @cols, 4; # first five columns only identify the variant

  # Prefix all column names with "am_"
  @cols = map { $_ =~ /^am_/ ? $_ : "am_" . $_ } @cols;
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

    # Prefix all column names with "am_"
    @cols = map { $_ =~ /^am_/ ? $_ : "am_" . $_ } @cols;
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
  die "\n  ERROR: No file specified\nTry using 'file=path/to/file.tsv.gz'\n"
     unless defined($file);
  $self->add_file($file);

  # Parse column names
  my $cols = $param_hash->{cols} || "am_pathogenicity:am_class";
  $self->_parse_colnames($cols);

  return $self;
}

sub feature_types {
  return [ 'Transcript' ];
}

sub get_header_info {
  my $self = shift;
  my %header;
  my @keys = @{ $self->{colnames} };

  my $suffix = "column from " . $self->{_files}[0];
  my @vals = map { $suffix } @keys;
  @header{ @keys } = @vals;

  # Custom headers
  $header{"am_pathogenicity"} = "AlphaMissense pathogenicity score; " . $suffix;
  $header{"am_class"} = "AlphaMissense pathogenicity prediction; " . $suffix;
  $header{"am_protein_variant"} = "Amino acid change used in AlphaMissense prediction; " . $suffix;
  $header{"am_uniprot_id"} = "Protein isoform used in AlphaMissense prediction; " . $suffix;
  $header{"am_transcript_id"} = "Transcript sequence in AlphaMissense prediction; " . $suffix;

  # Filter by user-selected columns
  %header = map { $_ => $header{$_} } @{ $self->{cols} };

  return \%header;
}

sub _aminoacid_changes_match {
  my ($self, $tva, $am_protein_var) = @_;

  my $am_aa_ref = substr($am_protein_var, 0, 1);
  my $am_aa_alt = substr($am_protein_var, -1);
  my $am_aa_pos = substr($am_protein_var, 1, -1);
  return 0 unless defined $am_aa_ref && defined $am_aa_alt && defined $am_aa_pos;

  my $vf_aa_ref = $tva->base_variation_feature_overlap->get_reference_TranscriptVariationAllele->peptide;
  my $vf_aa_alt = $tva->peptide;
  my $vf_aa_pos = $tva->base_variation_feature_overlap->translation_start;
  return 0 unless defined $vf_aa_ref && defined $vf_aa_alt && defined $vf_aa_pos;

  return $vf_aa_pos eq $am_aa_pos &&
         $vf_aa_ref eq $am_aa_ref &&
         $vf_aa_alt eq $am_aa_alt;
}

sub run {
  my ($self, $tva) = @_;

  # Only process missense variants
  return {} unless grep {$_->SO_term eq 'missense_variant'} @{$tva->get_all_OverlapConsequences};

  # Get allele
  my $allele = $tva->base_variation_feature->alt_alleles;

  my $vf = $tva->variation_feature;
  my @data = @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})};

  foreach (@data) {
    # Check if aminoacid changes match
    my $am_protein_var = $_->{result}->{am_protein_variant};
    next unless $self->_aminoacid_changes_match($tva, $am_protein_var);

    # Check if transcripts match
    if ($self->{transcript_match}) {
      my $am_transcript = $_->{result}->{am_transcript_id};
      my $am_transcript_no_suffix = $am_transcript;
      $am_transcript_no_suffix =~ s/\.[0-9]//g;
      next if $am_transcript_no_suffix ne $tva->transcript->{stable_id};
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
    return \%res if (@$matches);
  }
  return {};
}

sub parse_data {
  my ($self, $line) = @_;
  my ($chrom, $start, $ref, $alt, @vals) = split /\t/, $line;

  # VCF-like adjustment of mismatched substitutions for comparison with VEP
  if(length($alt) != length($ref)) {
    my $first_ref = substr($ref, 0, 1);
    my $first_alt = substr($alt, 0, 1);
    if ($first_ref eq $first_alt) {
      $start++;
      $ref = substr($ref, 1) || "-";
      $alt = substr($alt, 1) || "-";
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
