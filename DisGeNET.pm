=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

 DisGeNET

=head1 SYNOPSIS

 mv DisGeNET.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin DisGeNET,/path/to/disgenet/data.tsv.gz
  ./vep -i variations.vcf --plugin DisGeNET,/path/to/disgenet/data.tsv.gz,1

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds Variant-Disease-PMID associations from the DisGeNET database.
 It is available for GRCh38.

 Please cite the DisGeNET publication alongside the VEP if you use this resource:
 https://academic.oup.com/nar/article/48/D1/D845/5611674


 Running options:
 (Option 1) By default, this plugin includes three values = pmid:rsid:DisGeNET score.
 (Option 2) It can be run with the option '1' to also include diseases/phenotypes
 names reporting the Variant-pmid association.

 Output:
 The output includes three values: 
  - pmid of the publication reporting the Variant-Disease association
  - dbSNP variant Identifier
  - DisGENET score for the Variant-Disease association
  Example '20054638:rs10012:0.01'


 The following steps are necessary before running this plugin:
 This plugin uses file 'all_gene_disease_pmid_associations.tsv.gz'
 File can be downloaded from: https://www.disgenet.org/downloads

 gunzip all_gene_disease_pmid_associations.tsv.gz

 awk '($1 ~ /^snpId/ || $2 ~ /NA/) {next} {print $0 | "sort -k2,2 -k3,3n"}' 
 all_gene_disease_pmid_associations.tsv > all_variant_disease_pmid_associations_sorted.tsv

 bgzip all_variant_disease_pmid_associations_sorted.tsv
 tabix -s 2 -b 3 -e 3 all_variant_disease_pmid_associations_sorted.tsv.gz

 The plugin can then be run as default (Option 1):
 ./vep -i variations.vcf --plugin DisGeNET,/path/to/disgenet_file.tsv.gz

 or with an option to include disease/phenotype data (Option 2): 
 ./vep -i variations.vcf --plugin DisGeNET,/path/to/disgenet_file.tsv.gz,1

 Of notice: this plugin only matches the chromosome and the position in the
  chromosome, the alleles are not taken into account to append the DisGENET data.
  The rsid is provided in the output in order to help to filter the relevant data.


=cut
package DisGeNET;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();

  if(defined($self->params->[1])) {
    $self->{disease} = $self->params->[1];
  }

  return $self;
}

sub feature_types {
  return ['Feature','Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my %header;

  $header{'DisGeNET'} = 'Publications reporting the Variant-Disease association. Output includes three values pmid:rsid:score. pmid - is the pmid of the publication reporting the Variant-Disease association; rsid - dbSNP variant Identifier; score - DisGENET score for the Variant-Disease association';

  if($self->{disease}) {
    $header{'DisGeNET_disease'} = 'Name of the disease reporting the Variant-pmid association';
  }

  return \%header;
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;

  my $chr = $vf->{chr};

  my $end = $vf->{end};
  my $start = $vf->{start};
  ($start, $end) = ($end, $start) if $start > $end;

  my @data = @{$self->get_data($chr, $start, $end)};

  return {} unless(@data);

  my %hash;
  my @result;
  my %unique_values;
  my @diseases;
  my %unique_diseases;

  my $n_results = 1;

  foreach my $data_value (@data) {
    my $pmid = $data_value->{pmid};
    my $rsid = $data_value->{rsid};
    my $score = $data_value->{score};

    # Some publications are duplicated - same publications from different sources are in different rows
    # Check if pmid and rsid are not returned more than once
    if(!$unique_values{$pmid.':'.$rsid}++) {
      my $result_string = $pmid . ':' . $rsid . ":" . $score;
      push @result, $result_string;
    }

    if($self->{disease}) {
      my $disease_name = $data_value->{diseaseName};
      if(!$unique_diseases{$disease_name}++) {
        push @diseases, $disease_name;
      }
    }
  }

  $hash{'DisGeNET'} = join(',', @result);

  if($self->{disease}) {
    $hash{'DisGeNET_disease'} = join(',', @diseases);
  }

  return \%hash;
}

sub parse_data {
  my ($self, $line) = @_;

  # Data in file is: 
  # 'snpId, chromosome, position, DSI, DPI, diseaseId, diseaseName, diseaseType, diseaseClass,
  # diseaseSemanticType, score, EI, YearInitial, YearFinal, pmid, source'

  my @all_data = split /\t/, $line;

  return {
      rsid => $all_data[0],
      diseaseName => $all_data[6],
      score => $all_data[10],
      pmid => $all_data[14]
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
