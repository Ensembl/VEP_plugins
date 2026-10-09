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
 ./vep -i variations.vcf --plugin NMD

=head1 DESCRIPTION 
 This is a plugin for the Ensembl Variant Effect Predictor (VEP) 
 that predicts if a variant allows the 
 transcript escape nonsense-mediated mRNA decay based on certain rules.
 
 The rules are : 

 1. The variant location  falls in the last exon of the transcript.
                                        vvvv
      ES...EE..I.ES...EE.I.ES....EE.I.ES....EE 
 (ES= exon_start,EE = exon_end, I = intron, v = variant location)

 2. The variant location falls  50 bases upstream of the penultimate (second to the last ) exon.
                               vvv
      ES...EE..I.ES...EE.I.ES....EE.I.ES....EE 
 (ES= exon_start,EE = exon_end, I = intron, v = variant location)
     
 3. The variant falls in the first 100 coding bases in the transcript.
         vvv
      ..ES...EE..I.ES...EE.I.ES....EE.I.ES....EE 
 (ES= exon_start,EE = exon_end, I = intron, v = variant location)
  
 4. If the variant is in an intronless transcript, meaning only one exon exist in the transcript. 
 
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
  return {
    NMD => "Nonsense-mediated mRNA decay escaping variants prediction"};
}

sub run {
  my $self = shift; 
  my $tva = shift;
  # using stop_gained for now as advised. 
  return {} unless grep {$INCLUDE_SO{$_->SO_term}} @{$tva->get_all_OverlapConsequences};
  
  my $tr = $tva->transcript;
  my $tv = $tva->transcript_variation;
  # position of the variant in respect to the coding sequence. 
  my $variant_coding_region = $tv->cds_end if (defined $tv->cds_end);
  # checking for if the transcript is an intronless transcript 
  my $introns = $tr->get_all_Introns;
  my $exons = $tr->get_all_Exons;
  my $intron_count = scalar @{$introns};

  my $vf = $tva->variation_feature;
  
  # to check what position the transcript variation coding sequence end is at 
  # former method gave room for a bug 
  my $position_check = 1 if ( defined ($variant_coding_region) && $variant_coding_region <= 101) ;

  
 
  # to check if the variant is not on the coding sequence.
  if (!defined ($tv->cds_end) || !defined ($tv->cds_start)){
    return {};
  }


  # if statement to check if any of the rules is true
  if ( ( defined ($position_check) && $position_check == 1)
      || in_last_exon($vf, $exons)
      || in_penultimate_exon_50bp_end($vf, $exons, $tr->strand)
      || $intron_count == 0
     ) {
    return {
      NMD => $TERM{1},
    };
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

1;
