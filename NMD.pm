=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute
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
 that predicts if a stop_gained variant allows the 
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
 Identifying Genes Whose Mutant Transcripts Cause Dominant Disease Traits by Potential Gain-of-Function Alleles (Coban-Akdemir, 2018)
 The rules and impact of nonsense-mediated mRNA decay in human cancers (Lindeboom, 2016)


=cut 

package NMD;

use strict;
use warnings;

use base  qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

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
  return {} unless grep {$_->SO_term eq 'stop_gained'} @{$tva->get_all_OverlapConsequences};
  
  my $tr = $tva->transcript;
  my $tv = $tva->transcript_variation;
  # position of the variant in respect to the coding sequence. 
  my $variant_coding_region = $tv->cds_end;
  # checking for if the transcript is an intronless transcript 
  my @introns = $tr->get_all_Introns; 
  my $number = scalar @introns;
  # to check for if the variant location falls in the last exon
  my $check = $self->variant_exon_check($tva);

  # if statement to check if any of the rules is true
  if ($variant_coding_region < 101 || $check || $number == 0 ) {
    return {NMD => $TERM{1}};
  }
  else {
    return {};
  }
}

=head2 variant_exon_check 
  
  Args        : $tva (TranscriptionVariationAllele)
  Description : Checks if the variant location falls :
                  - in the last exon 
                  -  50bp towards the end of the second to the last exon(if it exists). 
                If it does, It returns 1 else it returns 0.
  Returntype  : Boolean
  Exceptions  : None
  Caller      : run 
  Status      : Stable

=cut
sub variant_exon_check {
  my $self = shift;
  my $tva = shift; 
  my $vf = $tva->variation_feature;
  my $tr = $tva->transcript;
  my $vf_end = $vf->seq_region_end;
  my @exons = @{ $tr->get_all_Exons }; 
  # to get the last exon 
  my $last_exon = $exons[-1]; 
  # to get the second to the last exon (penultimate exon)
  my $second_last_exon = $exons[-2];
  my $coding_region_end = $last_exon->coding_region_end($tr);
  my $coding_region_start = $last_exon->coding_region_start($tr);
  if (defined($coding_region_start) && $vf_end >= $coding_region_start && $vf_end <= $coding_region_end){
    return 1; 
  }
  # to check if the second to the last exon exists 
  if (defined($second_last_exon)){
    my $coding_region_end = $second_last_exon->coding_region_end($tr);
    if (defined($coding_region_end) && $vf_end >= $coding_region_end - (51) ){
      return 1; 
    }
  }
  
  return 0; 

}




1;

