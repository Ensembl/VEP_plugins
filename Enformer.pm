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
 Enformer

=head1 SYNOPSIS
 
 mv Enformer.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin Enformer,file=Enformer_grch38.vcf.gz


=head1 DESCRIPTION
 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that adds pre-calculated Enformer predictions of variant impact on chromatin and gene expression.

 The predictions have been aggregated across all 896 spatial bins to generate 5313 features corresponding to track prediction changes in differnet assays and cell types.

 This plugin is available for GRCh37 and GRCh38
 
 Please cite the Enformer publication alongside the VEP if you use this resource:
 https://www.nature.com/articles/s41592-021-01252-x
 
 GRCh38 scores were lifted over using CrossMap from the Enformer scores available here - https://console.cloud.google.com/storage/browser/dm-enformer/variant-scores/1000-genomes/enformer
 
 Enformer scores can be downloaded from https://ftp.ensembl.org/pub/current_variation/Enformer for GRCh37 and GRCh38.

 The plugin can then be run as default to retrieve SAD (SNP Activity Difference (SAD) and SAR (Same as SAD, by computing np.log2(1 + model(alternate_sequence)) - np.log2(1 + model(reference_sequence)) scores from Enforme :
 ./vep -i variations.vcf --assembly GRCh38 --plugin Enformer,file=/path/to/Enformer/data.vcf.gz

 or run with option to only retrieve the SAD (SNP Activity Difference (SAD) scores - main variant effect score computed as model(alternate_sequence) - model(reference_sequence) score 
 ./vep -i variations.vcf --assembly GRCh38 --plugin Enformer,file=/path/to/Enformer/data.vcf.gz,SAD=1 
 
 or run with option to only retrieve the SAR (Same as SAD, by computing np.log2(1 + model(alternate_sequence)) - np.log2(1 + model(reference_sequence)) score 
 ./vep -i variations.vcf --assembly GRCh38 --plugin Enformer,file=/path/to/Enformer/data.vcf.gz,SAR=1 

 or run with option to also retrieve the principal component scores which are a reduced representation of a much bigger vector of the SAD and SAR after using principal component analysis (PCA)
 ./vep -i variations.vcf --assembly GRCh38 --plugin Enformer,file=/path/to/Enformer/data.vcf.gz,PC=1 


 The tabix utility must be installed in your path to use this plugin.
 Check https://github.com/samtools/htslib.git for instructions.


=cut

package Enformer;

use strict;
use warnings; 

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();
  
  my $params = $self->params_to_hash();
  my $file;
   
  die "File needs to be specified to run the Enformer plugin. File can be downloaded from  https://ftp.ensembl.org/pub/current_variation/Enformer \n" if  (!%$params);
  if (!keys %$params) {
    $file = $self->params->[0];
    $params->{file} = $file;
  } else {
    $file = $params->{file};
    $self->{SAD} = $params->{SAD} if (defined ($params->{SAD}));
    $self->{SAR} = $params->{SAR} if (defined ($params->{SAR}));
    $self->{PC}  = $params->{PC} if (defined ($params->{PC}));
  } 

  delete $params->{SAD} if defined $params->{SAD} && $params->{SAD} eq '0';
  delete $params->{SAR} if defined $params->{SAR} && $params->{SAR} eq '0';
  delete $params->{PC} if defined $params->{PC} && $params->{PC} eq '0';
  
  $self->{params} = $params;
  $self->add_file($file);
  

  return $self;

}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my %header; 

  if (keys(%{$self->{params}}) == 1)  {
    $header{"Enformer_SAD"} = "Predictions of variant impact on gene expression. SNP Activity Difference (SAD) scores";
    $header{"Enformer_SAR"} = "Predictions of variant impact on gene expression. SNP Activity Difference logarithm computing";
  }

  if ($self->{SAD}) {
    $header{"Enformer_SAD"} = "Predictions of variant impact on gene expression. SNP Activity Difference (SAD) scores";
  }

  if ( $self->{SAR}) {
    $header{"Enformer_SAR"} = "Predictions of variant impact on gene expression. SNP Activity Difference logarithm computing";
  }
  
  if ($self->{PC}) {
    $header{"Enformer_SAD"} = "Predictions of variant impact on gene expression. SNP Activity Difference (SAD) scores";
    $header{"Enformer_SAR"} = "Predictions of variant impact on gene expression. SNP Activity Difference logarithm computing";
    $header{"Enformer_PC"} = "Predictions of variant impact on gene expression. Principal components of variant-effect scores";
  }

  return \%header;
}

sub run {
  my ($self, $tva) = @_;
  
  my $vf = $tva->variation_feature;
  my $allele = $tva->variation_feature_seq;

   # get allele
  my $alt_alleles = $tva->base_variation_feature->alt_alleles;
  my $ref_allele = $vf->ref_allele_string;
  
  my ($vf_start, $vf_end) = ($vf->{start}, $vf->{end});
  ($vf_start, $vf_end) = ($vf_end, $vf_start) if ($vf_start > $vf_end);
  
  my @data = @{
    $self->get_data(
      $vf->{chr},
      $vf_start,
      $vf_end 
    )
  };
  return {} unless(@data);
  
  # to account for insertions and deletions, need to do this search differently
  
  foreach my $variant (@data) {
    my $matches = get_matched_variant_alleles(
      {
        ref    => $ref_allele,
        alts   => $alt_alleles,
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
       ref  => $variant->{ref},
       alts => [$variant->{alt}],
       pos  => $variant->{start},
      }
    );
    return $variant->{result} if (@$matches) && keys(%{$self->{params}}) == 1;
    
    return $variant->{result_PC} if (@$matches) && $self->{PC};

    return {Enformer_SAD => $variant->{result}->{Enformer_SAD}} if (@$matches) && $self->{SAD};
    
    return {Enformer_SAR => $variant->{result}->{Enformer_SAR}} if (@$matches) && $self->{SAR};
    
  }

  return {};

}


sub parse_data {
  my ($self, $line) = @_;
  my ($chr_data, $start, $snp, $ref, $alt, $qual, $filter, $data) = split("\t", $line);
  
  my ($chr) = $chr_data =~ /chr(.+)/; # this is because the chromosome is chr1 etc, to retrieve just the 1

  my @data_splitted = split("SAD", $data); #splitting data in the vcf file based on the format using SAD
 
  my $POC = $data_splitted[0];
  
  my @other_scores = split(";", $data_splitted[1]); # to get SAD and SAR out of the scores 
  my $SAD = $other_scores[0] =~ s/=//r;  # SAD comes first 
  my $SAR = $other_scores[1] =~ s/.+=//r; # SAR is second 

  return {
    chr => $chr,
    start => $start,
    ref => $ref,
    alt => $alt,
    result_PC => {
      Enformer_SAD => $SAD,
      Enformer_SAR => $SAR,
      Enformer_PC  => $POC
    }, 
    result => {
      Enformer_SAD => $SAD,
      Enformer_SAR => $SAR
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