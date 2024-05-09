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

 AVADA

=head1 SYNOPSIS
./vep -i variations.vcf --plugin AVADA,file=path/to/file
./vep -i variations.vcf --plugin AVADA,file=path/to/file,feature_match_by=<gene_symbol|ensembl_gene_id|refseq_id>

=head1 DESCRIPTION

Automatic VAriant evidence DAtabase is a novel machine learning tool that uses natural language processing 
to automatically identify pathogenic genetic variant evidence in full-text primary literature about 
monogenic disease and convert it to genomic coordinates.

Please cite the AVADA publication alongside the VEP if you use this resource:
https://pubmed.ncbi.nlm.nih.gov/31467448/

NB: The plugin currently does not annotate for downstream_gene_variant and upstream_gene_variant.

Pre-requisites
1) AVADA data is available for GRCh37 and can be downloaded from: 
http://bejerano.stanford.edu/AVADA/avada_v1.00_2016.vcf.gz

wget http://bejerano.stanford.edu/AVADA/avada_v1.00_2016.vcf.gz

2) The file needs to be tabix indexed. You can do this by following commands:

gzip -d avada_v1.00_2016.vcf.gz
bgzip avada_v1.00_2016.vcf 
tabix avada_v1.00_2016.vcf.gz

3) As you have already noticed, tabix utility must be installed in your path to use this plugin. 

The plugin can then be run to retrieve AVADA annotations. 
By default, the variants are matched with the HGNC gene symbol
./vep -i variations.vcf --plugin AVADA,file=path/to/file

The output always includes one of the following columns depending on the option passed:
- `AVADA_PMID`: PubMed ID evidence for the variant as reported by AVADA
- `AVADA_PMID_WITH_VARIANT`: PubMed ID evidence for the variant as reported by AVADA along with the original variant string
- `AVADA_PMID_WITH_FEATURE`: PubMed ID evidence for the variant as reported by AVADA along with feature id
- `AVADA_PMID_WITH_FEATURE_AND_VARIANT`: PubMed ID evidence for the variant as reported by AVADA along with feature id and original variant string

The plugin can optionally be run by specifying the feature to match with.

In order to match by HGNC gene symbol:
./vep -i variations.vcf --plugin AVADA,file=path/to/file,feature_match_by=gene_symbol 

In order to match by Ensembl gene identifier :
./vep -i variations.vcf --plugin AVADA,file=path/to/file,feature_match_by=ensembl_gene_id

In order to match by RefSeq identifier :
./vep -i variations.vcf --plugin AVADA,file=path/to/file,feature_match_by=refseq_id

The plugin can also be run to report the original variant string reported in the publication. 
./vep -i variations.vcf --plugin AVADA,file=path/to/file,original_variant_string=1


=cut

package AVADA;

use strict;
use warnings;
use List::MoreUtils qw(uniq);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub get_header_info {
  my $self = shift;
  my %header;

  $header{"AVADA_PMID"} = "PubMed ID evidence for the variant as reported by AVADA" if not $self->{original_variant_string} and not $self->{feature_match_by} ;  
  $header{"AVADA_PMID_WITH_VARIANT"} = "PubMed ID evidence for the variant as reported by AVADA along with original variant string. The string is of format: PMID\%variant_string" if $self->{original_variant_string} and not $self->{feature_match_by}; 
  $header{"AVADA_PMID_WITH_FEATURE"} = "PubMed ID evidence for the variant as reported by AVADA along with feature id. The string is of format: PMID\%feature_id" if not $self->{original_variant_string} and $self->{feature_match_by}; 
  $header{"AVADA_PMID_WITH_FEATURE_AND_VARIANT"} = "PubMed ID evidence for the variant as reported by AVADA along with feature id and original variant string. The string is of format: PMID\%feature_id\%variant_string" if $self->{original_variant_string} and $self->{feature_match_by}; 

  return \%header
}

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();
  my $param_hash = $self->params_to_hash();
  my $file = $param_hash->{file};
  die "\n  ERROR: No file specified\nTry using 'AVADA,file=path/to/file.tsv.gz'\n" unless defined($file);
  $self->add_file($file);
  $self->{feature_match_by} = $param_hash->{feature_match_by}; 
  $self->{original_variant_string} = $param_hash->{original_variant_string}; 

  if ($self->{feature_match_by} and $self->{feature_match_by} eq "refseq_id" ) {
    die "\n ERROR: Matching by Refseq ID requires the option --refseq or --merged \n" unless (defined($self->{config}->{refseq}) || defined($self->{config}->{merged}));
  }
  return $self;
}

sub run {
  my ($self, $tva) = @_;

  return {} if grep {$_->SO_term eq 'downstream_gene_variant' || $_->SO_term eq 'upstream_gene_variant'} @{$tva->get_all_OverlapConsequences};

  my $vf = $tva->variation_feature;
  my $tv = $tva->transcript_variation;

  my ($vf_start, $vf_end) = ($vf->{start}, $vf->{end});
  ($vf_start, $vf_end) = ($vf_end, $vf_start) if ($vf_start > $vf_end);

  my $transcript = $tva->transcript;
  my $alt_allele = $tva->variation_feature_seq;
  my %output;
  my $refseq_transcript;
  my $refseq_protein;
  my $feature_key; 
  my @data;

  if (!defined($self->{feature_match_by})){
    @data = @{$self->get_data($vf->{chr}, $vf_start, $vf_end)}; 
  }
  elsif ( $self->{feature_match_by} eq "gene_symbol"){
    @data = grep {
    $_->{AVADA_GENE_SYMBOL} eq $transcript->{_gene_symbol}
    }@{$self->get_data($vf->{chr}, $vf_start, $vf_end)};
    $feature_key = "AVADA_GENE_SYMBOL";
  }
  elsif ( $self->{feature_match_by} eq "ensembl_gene_id" ){
    @data = grep {
    $_->{AVADA_ENSEMBL_ID} eq $transcript->{_gene_stable_id}
    }@{$self->get_data($vf->{chr}, $vf_start, $vf_end)};
    $feature_key = "AVADA_ENSEMBL_ID";
  }
  elsif ( $self->{feature_match_by} eq "refseq_id" ){
    $refseq_transcript = $transcript->{stable_id};
    $refseq_protein = $transcript->translation->{stable_id} if defined($transcript->translation);
    $refseq_protein =~ s/cds-// if $self->{config}->{database} == 1;
    @data = grep {
      split(".",$_->{AVADA_REFSEQ_ID}) eq split(".",$refseq_transcript) || 
      $_->{AVADA_REFSEQ_ID} eq $refseq_protein 
      } @data;
    $feature_key = "AVADA_REFSEQ_ID";
  }
  else{
    die("ERROR: feature_match_by can only take one of the options gene_symbol|ensembl_gene_id|refseq_id");
  }
  

  return {} unless scalar @data;

  my $pmid_string;
  my $output_key;
  my %seen;
  foreach my $data_value (uniq @data) {
    next unless $alt_allele eq $data_value->{AVADA_ALT};
    my $pmid_variant;
    my $pmid_variant_key;
    # Output (except json) is in string format (Eg: "PMID%feature_id%variant_string") 
    if (not $self->{config}->{output_format} eq 'json')
    {
      $pmid_variant = [$data_value->{AVADA_PMID}];
      push @$pmid_variant,$data_value->{$feature_key} if $self->{feature_match_by};
      push @$pmid_variant,$data_value->{AVADA_VARIANT_STRING} if $self->{original_variant_string};
      $pmid_variant = join('%', @$pmid_variant);
      $pmid_variant_key = $pmid_variant; 
    }
    else
    {
      $pmid_variant->{"avada_pmid"} = $data_value->{AVADA_PMID};
      $pmid_variant->{"avada_variant_string"} = $data_value->{AVADA_VARIANT_STRING} if $self->{original_variant_string};
      $pmid_variant->{"avada_feature_id"} = $data_value->{$feature_key} if $self->{feature_match_by};
      $pmid_variant_key =  join('%', values %{$pmid_variant});
    }
    next unless (! exists $seen{$pmid_variant_key});
    push @$pmid_string, $pmid_variant;
    $seen{$pmid_variant_key} = 1;

    $output_key = "AVADA_PMID_WITH_VARIANT" if $self->{original_variant_string} and not  $self->{feature_match_by};
    $output_key = "AVADA_PMID_WITH_FEATURE" if $self->{feature_match_by} and not $self->{original_variant_string};
    $output_key = "AVADA_PMID_WITH_FEATURE_AND_VARIANT" if $self->{feature_match_by} and $self->{original_variant_string};
    $output_key = "AVADA_PMID" if not $self->{feature_match_by} and not $self->{original_variant_string};; 
    
  }
  $output{$output_key} = $pmid_string if defined($output_key);
  return \%output;
}

sub parse_data {
  my ($self, $line) = @_;
  my ($chr, $start, $id, $ref, $alt, $x, $xx, $data) = split /\t/, $line;
  my ($pmid, $ensembl_gene_id, $gene_symbol, $refseq_id, $original_variant_string);
  my @data_split = split /;/, $data;
  foreach my $value (@data_split){
      $pmid = $value if $value =~ /PMID/;
      $ensembl_gene_id = $value if $value =~ /ENSEMBL_ID/; 
      $gene_symbol = $value if $value =~ /GENE_SYMBOL/; 
      $refseq_id = $value if $value =~ /REFSEQ_ID/; 
      $original_variant_string = $value if $value =~ /ORIGINAL_VARIANT_STRING/;
  }
  $pmid =~ s/PMID=//;
  $ensembl_gene_id =~ s/ENSEMBL_ID=//;
  $gene_symbol =~ s/GENE_SYMBOL=//;
  $refseq_id =~ s/REFSEQ_ID=//;
  $original_variant_string =~ s/ORIGINAL_VARIANT_STRING=//;
  return {
    AVADA_ALT => $alt,
    AVADA_PMID => $pmid,
    AVADA_ENSEMBL_ID => $ensembl_gene_id,
    AVADA_GENE_SYMBOL => $gene_symbol,
    AVADA_REFSEQ_ID => $refseq_id,
    AVADA_VARIANT_STRING => $original_variant_string
  };
}

1;
