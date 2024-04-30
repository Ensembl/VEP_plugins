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

 VARITY

=head1 SYNOPSIS

 mv VARITY.pm ~/.vep/Plugins
 ./vep -i variations.vcf --assembly GRCh37 --plugin VARITY,file=/path/to/varity_all_predictions.txt

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the pre-computed VARITY scores to predict pathogenicity of 
 rare missense variants to VEP output.

 Please cite the VARITY publication alongside the VEP if you use this resource:
 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8715197/

 Running options :

 VARITY scores can be downloaded using
 wget http://varity.varianteffect.org/downloads/varity_all_predictions.tar.gz

 The files can be tabix processed by :
 tar -xzvf varity_all_predictions.tar.gz 
 cat varity_all_predictions.txt | (head -n 1 && tail -n +2  | sort -t$'\t' -k 1,1 -k 2,2n) > varity_all_predictions_sorted.tsv
 sed '1s/.*/#&/'  varity_all_predictions_sorted.tsv > varity_all_predictions.tsv  # to add a # in the first line of the file 
 bgzip varity_all_predictions.tsv
 tabix -f -s 1 -b 2 -e 2 varity_all_predictions.tsv.gz


 Requirements:
 The tabix utility must be installed in your path to use this plugin.
 The --assembly flag is required to use this plugin.


=cut 


package VARITY;
 
use strict;
use warnings; 

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();
  
  my $params = $self->params_to_hash();

  my $file;
  if (!keys %$params) {
    $file = $self->params->[0];
    $params->{file} = $file;
  } else {
    $file = $params->{file};
  } 
  
  $self->add_file($file);
  
  my $assembly = $self->{config}->{assembly};

  if ($assembly ne "GRCh37") {
    die "Assembly is not GRCh37, VARITY only works with GRCh37. \n";
  }
  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return { VARITY_R => 'Pathogenicity score predicted by the VARITY_R model, specialized in rare missense variants (MAF < 0.5%)',
           VARITY_ER => 'Prediction by VARITY_ER model which was specialized in predicting extremely rare missense variants (No Allele count in gnomAD database)',
           VARITY_R_LOO => 'Pathogenicity score predicted by the VARITY_R model with the Leave-one-variant-out (LOO) strategy.',
           VARITY_ER_LOO => 'Pathogenicity score predicted by the VARITY_ER model with the Leave-one-variant-out (LOO) strategy.'
          };
}

sub run {
  my ($self, $tva) = @_;
  # only for missense variants
  return {} unless grep {$_->SO_term eq 'missense_variant'} @{$tva->get_all_OverlapConsequences};

  my $vf = $tva->variation_feature;

  my $allele = $tva->base_variation_feature->alt_alleles;


  my @data =  @{$self->get_data($vf->{chr}, $vf->{start}, $vf->{end})};
  return {} unless(@data);

  foreach my $variant (@data) {
    my $aa_ref = $tva->base_variation_feature_overlap->get_reference_TranscriptVariationAllele->peptide;
    my $aa_alt = $tva->peptide;
    my $aa_pos = $tva->base_variation_feature_overlap->translation_start;
    next unless defined $aa_pos && $aa_pos eq $variant->{aa_pos} &&
                defined $aa_ref && $aa_ref eq $variant->{aa_ref} &&
                defined $aa_alt && $aa_alt eq $variant->{aa_alt};

    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => $allele,
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
       ref  => $variant->{ref},
       alts => [$variant->{alt}],
       pos  => $variant->{start},
      }
    );
    #return $variant->{result} if ( (@$matches) && ($variant->{aa_alt} eq $tva->peptide) );
    return $variant->{result} if (@$matches) ;
  }

  return {};

}

sub parse_data {
  my ($self, $line) = @_;
  
  my @values = split /\t/, $line;

  my ($c, $pos, $ref, $alt, $p_vid, $aa_pos, $aa_ref, $aa_alt, $varity_r, $varity_er, $varity_r_loo, $varity_er_loo) = @values;
  

  return {
    chr => $c,
    start => $pos, 
    ref => $ref, 
    alt => $alt, 
    p_vid => $p_vid,
    aa_pos => $aa_pos, 
    aa_ref => $aa_ref, 
    aa_alt => $aa_alt, 
    result => {
      VARITY_R => $varity_r,
      VARITY_ER => $varity_er, 
      VARITY_R_LOO => $varity_r_loo, 
      VARITY_ER_LOO => $varity_er_loo
    }
  };

}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;

