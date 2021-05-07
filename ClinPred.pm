=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

     CLINPRED

=head1 SYNOPSIS
 mv Clinpred.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin CLINPRED,/path/to/Clinpred/data.txt.gz
=head1 DESCRIPTION
 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the CLINPRED(pathogenicity) score for variants to VEP output.

 Please cite the CLINPRED publication alongside the VEP if you use this resource:

CLINPRED SCORES CAN BE DOWNLOADED FROM 
https://sites.google.com/site/clinpred/products-services
 and can be tabix-processed by:
 
 for GRCh37:
 tabix -f -s 1 -b 2 -e 2 ClinPred.txt.gz 

 The tabix utility must be installed in your path to use this plugin.



=cut
package CLINPRED;
 

use strict;
use warnings; 



use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp)
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();

  return $self;
}

sub feature_types {
  return ['Transcript'];

}

sub get_header_info {
     # a bit confused about what name to call this since, it is not an abbreviation 
     return { ClinPred => "AN ENSEMBLE MACHINE LEARNING PREDICTION TOOL" };
}

sub run{
     my($self, $tva) = @_;

     # since it is only for missense variant, we do an if or unless statement 
     if (grep {$_->eq "missense_variant"}){
          return @{$tva->get_all_OverlapConsequences()};
     }

     else {
          return {}
     }

     my $vf = $tva->variation_feature();
     # get allele here from the vcf file 
     my $allele = $tva->variation_feature_seq();
}


sub parse_data {
     # we are parsing data using the header of the ClinPred.txt.gz file 
     # parse data is a method called from the base class BaseVepTabixPlugin

     my ($self, $line) = @_;
     # splitting the header to create a list of variables. 
     my ($chr, $start, $ref, $alt, $ClinPred_score ) = split("\t", $line);

     return {
          start => $start,
          end => $start,
          alt => $alt,
          result => {
               CLINPRED => $ClinPred_score,
          }
           };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}