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

=head1 DESCRIPTION
 This is a plugin for the Ensembl Variant Effect Predictor (VEP)
 that predicts if a variant allows the resulting isoform to escape
 nonsense-mediated mRNA decay, and categorizes those that do
 based on certain rules.

 Options are passed to the plugin as key=value pairs:
   report_category : (optional) boolean to enable/disable NMD category reporting. 1 (enabled) or 0 (disabled) (default: 0)

 The NMD escape categorisation rules are as listed below. Start-proximal rules (rule 3) apply to all variants.
 Rules based on the last coding exon-exon junction (rules 1, 2 and 4) do not apply to variants in transcripts with biotype `nonsense_mediated_decay`,
 as these transcripts are known to already have a PTC >= 50nt before the last coding exon-exon junction, thus missing the window of escape.

 1. The variant location falls in the last exon in the coding sequence of the transcript.
                                        vvvv
      ES...EE..I.ES...EE.I.ES....EE.I.ES....EE
    Category: `cds_E-1`
 (ES= exon_start,EE = exon_end, I = intron, v = variant location)

 2. The variant location falls within the last 50bp (3' end)
    of the penultimate exon in the coding sequence of the transcript.
                               vvv
      ES...EE..I.ES...EE.I.ES....EE.I.ES....EE
    Category: `cds_E-2_end-50`
 (ES= exon_start,EE = exon_end, I = intron, v = variant location)

 3. The variant falls in the first 100 coding bases in the transcript.
         vvv
      ..ES...EE..I.ES...EE.I.ES....EE.I.ES....EE
    Category: `cds_start+100`
 (ES= exon_start,EE = exon_end, I = intron, v = variant location)

 4. The variant falls anywhere in the coding sequence of a transcript with no introns in its coding sequence,
    i.e. the coding sequence is completely encapsulated in a single exon of the transcript (introns can exist in the UTRs).
    Category: `cds_intronless`

 The additional term NMD-escaping variant (nonsense-mediated mRNA decay escaping variants) will be added if the variant matches any of the rules.

 REFERENCES :
 - Identifying Genes Whose Mutant Transcripts Cause Dominant Disease Traits by Potential Gain-of-Function Alleles (Coban-Akdemir, 2018)
 - The rules and impact of nonsense-mediated mRNA decay in human cancers (Lindeboom, 2016)


=cut

package NMD;

use strict;
use warnings;

use base  qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my %INCLUDE_SO = map{$_ => 1} qw(stop_gained frameshift_variant splice_donor_variant splice_acceptor_variant);
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
    $header_info->{NMD_category} = "Nonsense-mediated mRNA decay escape category. One of:"
      ." 'E-1' (last exon), 'E-2_end-50' (last 50bp in the penultimate exon)"
      .", 'cds_start+100' (first 100bp of the coding sequence), 'intronless' (intronless transcript).";
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

  # Get the exons in the coding sequence region only.
  # UTR introns are irrelevant to determine NMD escape
  my $exons = $tr->get_all_translateable_Exons;


  # Skip NMD annotation if variant does not fall in the coding sequence
  if (!defined ($tv->cds_end) || !defined ($tv->cds_start)){
    return {};
  }

  my $tr_biotype = $tr->biotype;

  # Categorize the variant to NMD-escape classes
  my $nmd_escape_category;
  if( $tr_biotype ne 'nonsense_mediated_decay' && intronless_transcript_cds($exons) ){
    $nmd_escape_category = 'cds_intronless';
  }
  elsif( $tr_biotype ne 'nonsense_mediated_decay' && in_last_exon($vf, $exons) ){
    $nmd_escape_category = 'cds_E-1';
  }
  elsif( $tr_biotype ne 'nonsense_mediated_decay' && in_penultimate_exon_50bp_end($vf, $exons, $tr->strand) ){
    $nmd_escape_category = 'cds_E-2_end-50';
  }
  elsif( in_100bp_coding_start($tv) ){
    $nmd_escape_category = 'cds_start+100';
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

  my $last_50bp_start;
  my $last_50bp_end;
  # need to account for strand as if it is in reverse strand it should be reading from the start
  # and adding 51 to find the variant 50 bases upstream in penultimate exon
  if ($strand == -1){
    my $coding_region_start = $second_last_exon->start;
    # skip if coding region start could not be found
    if (!defined($coding_region_start)) {
      return 0;
    }

    $last_50bp_start = $coding_region_start;
    $last_50bp_end = $last_50bp_start + 51;
  }
  else {
    my $coding_region_end = $second_last_exon->end;
    # skip if coding region end could not be found
    if (!defined($coding_region_end)) {
      return 0;
    }

    $last_50bp_end = $coding_region_end;
    $last_50bp_start = $last_50bp_end - 51;
  }

  my $vf_end = $vf->seq_region_end;

  return $last_50bp_start <= $vf_end && $vf_end <= $last_50bp_end;
}

=head2 in_100bp_coding_start

  Args        : $tv (TranscriptVariation)
  Description : Checks if the variant location falls in the first 100bp (5' end)
                of the transcript's coding sequence.
                If so, returns 1 else it returns 0.
  Returntype  : Boolean
  Exceptions  : None
  Caller      : run
  Status      : Development

=cut
sub in_100bp_coding_start {
  my $tv = shift;

  # position of the variant in respect to the coding sequence.
  my $variant_coding_region = $tv->cds_end // undef;

  if(!defined($variant_coding_region)) {
    return 0;
  }

  return $variant_coding_region <= 101;
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
