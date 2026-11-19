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
 NMD

=head1 SYNOPSIS

  mv NMD.pm  ~/.vep/Plugins

  # Default usage
  ./vep -i variations.vcf --plugin NMD
  # Enable NMD category reporting (as NMD_category)
  ./vep -i variations.vcf --plugin NMD,report_category=1
  # Enable NMD reporting in the extended start region (in addition to the defaults)
  ./vep -i variations.vcf --plugin NMD,cds_start_boundary=200
  # Enable all optional NMD category rules
  ./vep -i variations.vcf --plugin NMD,cds_start_boundary=200,apply_categories=all
  # Select specific categories to report (quote the argument to protect the pipe from the shell)
  ./vep -i variations.vcf --plugin 'NMD,apply_categories=cds_intronless|cds_long_exon'

=head1 DESCRIPTION
 This is a plugin for the Ensembl Variant Effect Predictor (VEP)
 that predicts if a variant allows the resulting isoform to escape
 nonsense-mediated mRNA decay, and categorizes those that do
 based on certain rules.

 Options are passed to the plugin as key=value pairs:
   report_category :    (optional) boolean to enable/disable NMD category reporting. 1 (enabled) or 0 (disabled) (default: 0)
   cds_start_boundary : (optional) integer defining the boundary to which variants are considered to match start-proximal NMD escape rules (`cds_start+100` and `cds_start_extended`)
                        Must be between 100 and 200 (inclusive). When set > 100, apply_categories is implied to include `cds_start_extended`.
                        This will add cds_start_extended to the default ruleset if apply_categories is not set and error if 'apply_categories' is set but cds_start_extended is not included (default: 100)
   apply_categories :   (optional) pipe-separated list of categories to apply, or 'all' to apply all categories.
                        Enabling cds_start_extended requires cds_start_boundary to be set > 100.
                        (default: cds_intronless|cds_E-1|cds_E-2_end-50|cds_start+100)
                        See the rules below for the available category names.

 The NMD escape categorisation rules are as listed below. However, rules apply differently depending on the transcript biotype.
 * Start-proximal rules (rules 4-5) apply to all variants.
 * Rules based on the last coding exon-exon junction (rules 1-3) do not apply to variants in transcripts with biotype `nonsense_mediated_decay`,
   as these transcripts are known to already have a PTC >= 50nt before the last coding exon-exon junction, thus missing the window of escape.
 * Rule 6 (long exon) applies to all variants, regardless of transcript biotype.
   However, in transcripts with biotype `nonsense_mediated_decay`, the first downstream exon-exon junction past the CDS end
   can affect the escape potential of the variant, as 3' UTR introns in those transcripts are NMD-sensitive (unlike 3' UTR introns in other transcripts).
   As a such, the first downstream exon-exon junction past the CDS end
   is included in the evaluation of the rule for `nonsense_mediated_decay` transcripts.
 Rules are evaluated in documented order and only a single category will be reported (i.e. variants that match multiple rules
 will be annotated with the first matching category).

 Rules:
 1. The variant falls anywhere in the coding sequence of a transcript with no introns in its coding sequence,
    i.e. the coding sequence is completely encapsulated in a single exon of the transcript (introns can exist in the UTRs).
    Category: `cds_intronless`

 2. The variant location falls in the last exon in the coding sequence of the transcript.
                                        vvvv
      ES...EE..I.ES...EE.I.ES....EE.I.ES....EE
    Category: `cds_E-1`
 (ES= exon_start,EE = exon_end, I = intron, v = variant location)

 3. The variant location falls within the last 50bp (3' end)
    of the penultimate exon in the coding sequence of the transcript.
                               vvv
      ES...EE..I.ES...EE.I.ES....EE.I.ES....EE
    Category: `cds_E-2_end-50`
 (ES= exon_start,EE = exon_end, I = intron, v = variant location)

 4. The variant falls in the first 100 coding bases in the transcript.
         vvv
      ..ES...EE..I.ES...EE.I.ES....EE.I.ES....EE
    Category: `cds_start+100`
 (ES= exon_start,EE = exon_end, I = intron, v = variant location)

 5. The variant falls between coding bases 101 and `cds_start_boundary` (max 200) in the transcript.
                     vvv
      ..ES...EE..I.ES...EE.I.ES....EE.I.ES....EE
    This causes gradual low-efficiency NMD, with NMD efficiency increasing
    as the distance (in CDS bases) from the CDS start increases.
    This gives variants in this category a potential for partial NMD escape.
    Category: `cds_start_extended`
 (ES= exon_start,EE = exon_end, I = intron, v = variant location)

 6. The variant falls in a long exon within the CDS, with at least 200 bases
    between the variant and the next downstream in-CDS exon-exon junction.
    For NMD transcripts, the first downstream exon-exon junction past the CDS end
    is included in the evaluation of this rule.
    This causes gradual low-efficiency NMD, with NMD efficiency decreasing
    as the distance to the next downstream exon-exon junction increases.
    This gives variants in this category a potential for partial NMD escape.
    Category: `cds_long_exon`

 The additional term NMD-escaping variant (nonsense-mediated mRNA decay escaping variants) will be added if the variant matches any of the rules.

 Note: all rules are currently implemented based on the location of the variant, not the location of the premature stop codon.
 While these are sufficiently close in most cases, this is not a perfect implementation for large indels (a known shortcoming).
 This will not be fixed in the current version of the plugin, but will be reimplemented in the upcoming python reimplementation of ensembl-VEP.
 For more details on the ensembl-vep python reimplementation, see: https://github.com/Ensembl/ensembl-vep-py/.

 REFERENCES :
 - Identifying Genes Whose Mutant Transcripts Cause Dominant Disease Traits by Potential Gain-of-Function Alleles (Coban-Akdemir, 2018)
 - The rules and impact of nonsense-mediated mRNA decay in human cancers (Lindeboom, 2016)
 - Widespread 3′ UTR splicing regulates expression of oncogene transcripts through multiple mechanisms (Riley, 2025)


=cut

package NMD;

use strict;
use warnings;

use base  qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my %INCLUDE_SO = map{$_ => 1} qw(stop_gained);
# NMD escaping categories, in order of priority (matching rules order as described above)
my @CATEGORIES = ('cds_intronless', 'cds_E-1', 'cds_E-2_end-50',
                  'cds_start+100', 'cds_start_extended', 'cds_long_exon');
my %VALID_CATEGORY = map { $_ => 1 } @CATEGORIES;
my %TERM = (
  1 => "NMD_escaping_variant"
);

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  my $header_info = {
    NMD => "Nonsense-mediated mRNA decay escaping variants prediction"
  };

  if($self->{report_category}) {
    $header_info->{NMD_category} = "Nonsense-mediated mRNA decay escape category. One of:";

    if($self->{selected_categories}->{'cds_E-1'}) {
      $header_info->{NMD_category} .= " 'cds_E-1' (last exon),";
    }
    if($self->{selected_categories}->{'cds_E-2_end-50'}) {
      $header_info->{NMD_category} .= " 'cds_E-2_end-50' (last 50bp in the penultimate exon),";
    }
    if($self->{selected_categories}->{'cds_start+100'}) {
      $header_info->{NMD_category} .= " 'cds_start+100' (first 100bp of the coding sequence),";
    }
    if($self->{selected_categories}->{'cds_start_extended'}) {
      $header_info->{NMD_category} .= " 'cds_start_extended' (between bases 101 and ".$self->{cds_start_boundary}." of the coding sequence),";
    }
    if($self->{selected_categories}->{'cds_long_exon'}) {
      $header_info->{NMD_category} .= " 'cds_long_exon' (at least 200 bases upstream of an NMD-active EJC-bearing exon-exon junction),";
    }
    if($self->{selected_categories}->{'cds_intronless'}) {
      $header_info->{NMD_category} .= " 'cds_intronless' (intronless coding sequence),";
    }

    chop($header_info->{NMD_category});
  }

  return $header_info;
}

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $param_hash = $self->params_to_hash();

  my $report_category = $param_hash->{report_category};

  if ($report_category) {
    $self->{report_category} = 1;
  }
  else {
    $self->{report_category} = 0;
  }

  my $cds_start_boundary_str = $param_hash->{cds_start_boundary} // "100";
  my $cds_start_boundary = int($cds_start_boundary_str);

  if ( !$cds_start_boundary || "".$cds_start_boundary ne "".$cds_start_boundary_str ) {
    die "ERROR: cds_start_boundary is not an integer. Please provide an integer between 100 and 200 (inclusive).\n";
  }

  if ( $cds_start_boundary < 100 || 200 < $cds_start_boundary ) {
    die "ERROR: cds_start_boundary not between 100 and 200 (inclusive). Please provide an integer between 100 and 200 (inclusive).\n";
  }

  $self->{cds_start_boundary} = $cds_start_boundary;

  my $category_param = $param_hash->{apply_categories} // 'cds_intronless|cds_E-1|cds_E-2_end-50|cds_start+100';
  chomp($category_param);

  my @raw_categories = $category_param eq 'all' ? @CATEGORIES
                 : split(/\|/, $category_param, -1);
  my @invalid_categories = grep { !$VALID_CATEGORY{$_} } @raw_categories;
  if (!@raw_categories || @invalid_categories) {
    die "ERROR: Invalid apply_categories value. Use 'all' or a pipe-separated list of: "
      . join(', ', @CATEGORIES) . ".\n";
  }
  $self->{selected_categories} = { map { $_ => 1 } @raw_categories };

  ### Cross-check provided options

  # if cds_start_boundary > 100, cds_start_extended should be enabled
  if ($self->{cds_start_boundary} > 100) {
    # Enable cds_start_extended if apply_categories was not defined by user
    if( !$param_hash->{apply_categories} ) {
      $self->{selected_categories}{'cds_start_extended'} = 1;
    }
    else{
      # Error if user explicitly disabled cds_start_extended (providing incomaptible input args)
      if( !$self->{selected_categories}{'cds_start_extended'} ) {
        die "ERROR: cds_start_boundary > 100 implies cds_start_extended to be enabled."
          ." Please enable 'cds_start_extended' through the apply_categories option or remove the cds_start_boundary option.\n";
      }
    }
  }
  # Otherwise, cds_start_extended should be disabled
  else{
    if( $self->{selected_categories}{'cds_start_extended'} ) {
      die "ERROR: cds_start_boundary <= 100 implies cds_start_extended to be disabled."
        ." Please disable 'cds_start_extended' through the apply_categories option or specify the cds_start_boundary option to be > 100.\n";
    }
  }


  return $self;
}

sub run {
  my $self = shift;
  my $tva = shift;

  # Skip NMD annotation if transcript variant consequence is not subject to NMD
  return {} unless grep {$INCLUDE_SO{$_->SO_term}} @{$tva->get_all_OverlapConsequences};

  my $tr = $tva->transcript;
  my $tv = $tva->transcript_variation;
  my $vf = $tva->variation_feature;

  # Get the exons of the transcript
  my $cds_exons = $tr->get_all_translateable_Exons;  # in the coding sequence region only
  my $all_exons = $tr->get_all_Exons;  # in 5' UTR, coding sequence region and 3' UTR

  my $cds_and_downstream_exons = [];  # in the coding sequence region and 3' UTR

  # Filter all exons down to the ones including the CDS start and beyond (CDS + 3' UTR).
  # This assumes both @$cds_exons and @$all_exons are ordered 5' to 3'.
  if (@{$cds_exons}) {
    my $first_cds_exon = $cds_exons->[0];
    my $found_cds = 0;

    $cds_and_downstream_exons = [grep {
      $found_cds ||= (
        $_->start <= $first_cds_exon->end &&
        $_->end   >= $first_cds_exon->start
      );
    } @{$all_exons}];
  }

  # Select the exons to use for the long-exon rule
  my $long_exon_exons;
  if ($tr->biotype eq 'nonsense_mediated_decay') {
    $long_exon_exons = $cds_and_downstream_exons;
  } else {
    $long_exon_exons = $cds_exons;
  }


  # Skip NMD annotation if variant does not fall in the coding sequence
  if (!defined ($tv->cds_end) || !defined ($tv->cds_start)){
    return {};
  }

  my $tr_biotype = $tr->biotype;

  # Categorize the variant to NMD-escape classes
  my $nmd_escape_category;
  if( $self->{selected_categories}->{'cds_intronless'} && $tr_biotype ne 'nonsense_mediated_decay' && intronless_transcript_cds($cds_exons) ){
    $nmd_escape_category = 'cds_intronless';
  }
  elsif( $self->{selected_categories}->{'cds_E-1'} && $tr_biotype ne 'nonsense_mediated_decay' && in_last_exon($vf, $cds_exons) ){
    $nmd_escape_category = 'cds_E-1';
  }
  elsif( $self->{selected_categories}->{'cds_E-2_end-50'} && $tr_biotype ne 'nonsense_mediated_decay' && in_penultimate_exon_50bp_end($vf, $cds_exons, $tr->strand) ){
    $nmd_escape_category = 'cds_E-2_end-50';
  }
  elsif( $self->{selected_categories}->{'cds_start+100'} && in_coding_start($tv, 100) ){
    $nmd_escape_category = 'cds_start+100';
  }
  elsif( $self->{selected_categories}->{'cds_start_extended'} && !in_coding_start($tv, 100) && in_coding_start($tv, $self->{cds_start_boundary}) ){
    $nmd_escape_category = 'cds_start_extended';
  }
  elsif( $self->{selected_categories}->{'cds_long_exon'} && in_cds_long_exon($vf, $long_exon_exons, $tr->strand) ){
    $nmd_escape_category = 'cds_long_exon';
  }


  # Return NMD annotation if variant matches any of the rules
  if ( defined($nmd_escape_category) ) {
    my $annot = {NMD => $TERM{1}};
    if( $self->{report_category} ) {
      $annot->{NMD_category} = $nmd_escape_category;
    }
    return $annot;
  }
  else {
    return {};
  }


}

=head2 in_cds_long_exon

  Args        : ($vf (VariationFeature), \@exons (array of exons), $strand)
  Description : Checks whether there is at least 200 bases between the variant
                and the next downstream in-CDS exon-exon junction in the transcript.
  Returntype  : Boolean
  Caller      : run

=cut
sub in_cds_long_exon {
  my ($vf, $exons, $strand) = @_;

  my $position = $strand == 1 ? $vf->seq_region_end : $vf->seq_region_start;

  # Exons are ordered 5' to 3';
  # Last exon is considered not to have a following exon-exon junction (end of transcript)
  # or the following exon-exon junction is considered irrelevant (end of CDS),
  # so needs excluding in both cases.
  for (my $i = 0; $i < (scalar @{$exons}) - 1; $i++) {
    my $start = $exons->[$i]->start;
    my $end = $exons->[$i]->end;

    # skip if start or end could not be found
    if( !defined($start) || !defined($end) ) {
      next;
    }

    # skip if the variant does not fall in this exon
    if( !($start <= $position && $position <= $end) ){
      next;
    }

    return $strand == 1 ? $end - $position >= 200
                        : $position - $start >= 200;
  }

  return 0;
}

=head2 in_last_exon

  Args        : ($vf (VariationFeature), \@exons (array of Exons))
  Description : Checks if the variant location falls in the last exon
                If it does, It returns 1 else it returns 0.
  Returntype  : Boolean
  Exceptions  : None
  Caller      : run
  Status      : Development

=cut
sub in_last_exon {
  my $vf = shift;
  my $exons = shift;

  # Get the end of the variant
  my $vf_end = $vf->seq_region_end;

  # Get the last exon
  my $last_exon = $exons->[-1];
  my $coding_region_end = $last_exon->end;
  my $coding_region_start = $last_exon->start;

  # Skip if exon start or end could not be found
  if (!defined($coding_region_start) || !defined($coding_region_end)) {
    return 0;
  }

  # Check if the variant is in the last exon
  return $coding_region_start <= $vf_end && $vf_end <= $coding_region_end;
}

=head2 in_penultimate_exon_50bp_end

  Args        : ($vf (VariationFeature), \@exons (array of Exons))
  Description : Checks if the variant location falls in the 50bp towards the 3' end of the penultimate exon (if it exists).
                If so, returns 1 else it returns 0.
  Returntype  : Boolean
  Exceptions  : None
  Caller      : run
  Status      : Development

=cut
sub in_penultimate_exon_50bp_end {
  my $vf = shift;
  my $exons = shift;
  my $strand = shift;

  my $second_last_exon = $exons->[-2];

  if(!defined($second_last_exon)) {
    return 0;
  }

  my $exon_start = $second_last_exon->start;
  my $exon_end = $second_last_exon->end;
  # skip if coding region start or end could not be found
  if (!defined($exon_start) || !defined($exon_end)) {
    return 0;
  }

  my $last_50bp_start;
  my $last_50bp_end;
  # need to account for strand as if it is in reverse strand it should be reading from the start
  # and adding 51 to find the variant 50 bases upstream in penultimate exon
  if ($strand == -1){
    # 1-based inclusive positions
    $last_50bp_start = $exon_start;
    $last_50bp_end = $last_50bp_start + 49;
  }
  else {
    # 1-based inclusive positions
    $last_50bp_end = $exon_end;
    $last_50bp_start = $last_50bp_end - 49;
  }

  my $vf_end = $vf->seq_region_end;

  return $exon_start <= $vf_end && $vf_end <= $exon_end
         && $last_50bp_start <= $vf_end && $vf_end <= $last_50bp_end;
}

=head2 in_coding_start

  Args        : ($tv (TranscriptVariation), $limit (the coding start boundary (aka end)))
  Description : Checks if the variant location falls in the first $limit bp (5' end)
                of the transcript's coding sequence.
                If so, returns 1 else it returns 0.
  Returntype  : Boolean
  Exceptions  : None
  Caller      : run
  Status      : Development

=cut
sub in_coding_start {
  my $tv = shift;
  my $limit = shift;

  my $limit_int = int($limit);

  # Skip if limit is not defined or 0
  if (!defined($limit) || $limit == 0 || $limit_int != $limit) {
    if ($limit_int != $limit) {
      warn("in_coding_start: limit param must be an integer, evaluating as if limit is not defined.\n");
    }
    return 0;
  }

  # 1-based (inclusive) position of the variant relative to the coding sequence.
  my $variant_coding_region = $tv->cds_end // undef;

  if(!defined($variant_coding_region)) {
    return 0;
  }

  return $variant_coding_region <= $limit_int;
}

=head2 intronless_transcript_cds

  Args        : $tr (Transcript)
  Description : Checks if the transcript's coding sequence is intronless (single-exon).
                If so, returns 1 else it returns 0.
  Returntype  : Boolean
  Exceptions  : None
  Caller      : run
  Status      : Development

=cut
sub intronless_transcript_cds {
  my $exons = shift;

  my $exon_count = scalar @{$exons};

  return $exon_count == 1;
}

1;
