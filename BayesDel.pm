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

 BayesDel

=head1 SYNOPSIS

 mv BayesDel.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin BayesDel,file=/path/to/BayesDel/BayesDel_170824_addAF_all_scores.txt.gz

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the BayesDel scores to VEP output.

 BayesDel is a deleteriousness meta-score combining multiple deleteriousness predictors to create an overall score. It works for coding and non-coding variants,
 single nucleotide variants and small insertion/deletions.
 The range of the score is from -1.29334 to 0.75731.
 The higher the score, the more likely the variant is pathogenic.
 For more information please visit: https://fenglab.chpc.utah.edu/BayesDel/BayesDel.html

 Please cite the BayesDel publication alongside the Ensembl VEP if you use this resource:
 https://onlinelibrary.wiley.com/doi/full/10.1002/humu.23158


 BayesDel pre-computed scores can be downloaded from
 https://drive.google.com/drive/folders/1K4LI6ZSsUGBhHoChUtegC8bgCt7hbQlA
 Note: These files only contain pre-computed BayesDel scores for missense variants for assembly GRCh37.


 For GRCh37:
 > tar zxvf BayesDel_170824_addAF.tgz
 > rm *.gz.tbi
 > gunzip *.gz
 > for f in BayesDel_170824_addAF_chr*; do grep -v "^#" $f >> BayesDel_170824_addAF.txt; done
 > cat BayesDel_170824_addAF.txt | sort -k1,1 -k2,2n > BayesDel_170824_addAF_sorted.txt
 > grep "^#" BayesDel_170824_addAF_chr1 > BayesDel_170824_addAF_all_scores.txt
 > cat BayesDel_170824_addAF_sorted.txt >> BayesDel_170824_addAF_all_scores.txt
 > bgzip BayesDel_170824_addAF_all_scores.txt
 > tabix -s 1 -b 2 -e 2 BayesDel_170824_addAF_all_scores.txt.gz

 For GRCh38:
 Remap GRCh37 file


 The tabix utility must be installed in your path to use this plugin.

=cut

package BayesDel;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  my $param_hash = $self->params_to_hash();

  # Test tabix
  die "ERROR: tabix does not seem to be in your path\n" unless `which tabix 2>&1` =~ /tabix$/;

  # Get file
  if(defined($param_hash->{file})) {
    die "ERROR: BayesDel file not not found!\n" unless -e $param_hash->{file};
    $self->add_file($param_hash->{file});
  }

  # Get headers from file
  my $headers;
  open HEAD,"tabix -H $param_hash->{file} 2>&1 | ";
  while(<HEAD>) {
    chomp;
    $_ =~ s/^\#//;
    $headers = [split];
  }
  close HEAD;

  die "ERROR: Could not read headers from $param_hash->{file}\n" unless defined($headers) && scalar @{$headers};

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  return { BayesDel => 'Deleteriousness meta-score'};
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  my $transcript = $tva->transcript;

  my $end = $vf->{end};
  my $start = $vf->{start};
  ($start, $end) = ($end, $start) if $start > $end;

  my $ref_allele = $vf->ref_allele_string;
  my $alt_allele = $tva->variation_feature_seq;
  if ($vf->{strand} == -1) {
    reverse_comp(\$ref_allele);
    reverse_comp(\$alt_allele);
  }

  my ($res) = grep {
    $_->{start} == $start &&
    $_->{ref} eq $ref_allele &&
    $_->{alt} eq $alt_allele
  } @{$self->get_data($vf->{chr}, $start, $end)};

  return defined $res->{result} ? $res->{result} : {};
}

sub parse_data {
  my ($self, $line) = @_;
  my ($chr, $pos, $ref, $alt, $score) = split /\t/, $line;

  return {
    start => $pos,
    ref => $ref,
    alt => $alt,
    result => {
      BayesDel => $score
    }
  };
}

sub get_start {
  return $_[1]->{start};
}

1;
