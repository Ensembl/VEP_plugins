=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

 Please email comments or questions to the public Ensembl
 developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl helpdesk at
 <https://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

pLI - Add pLI score to the VEP output 

=head1 SYNOPSIS

  mv pLI.pm ~/.vep/Plugins
  mv pLI_values.txt ~/.vep/Plugins
  ./vep -i variants.vcf --plugin pLI

=head1 DESCRIPTION


  A VEP plugin that adds the probabililty of a gene being 
  loss-of-function intolerant (pLI) to the VEP output.
  
  Lek et al. (2016) estimated pLI using the expectation-maximization 
  (EM) algorithm and data from 60,706 individuals from 
  ExAC (http://exac.broadinstitute.org/about). The closer pLI is to 1, 
  the more likely the gene is loss-of-function (LoF) intolerant. 
  
  Note: the pLI was calculated using a representative transcript and
  is reported by gene in the plugin.

  The data for the plugin is provided by Kaitlin Samocha and Daniel MacArthur. 
  See https://www.ncbi.nlm.nih.gov/pubmed/27535533 for a description 
  of the dataset and analysis.

  The pLI_values.txt file is found alongside the plugin in the 
  VEP_plugins GitHub repository. The file contains the fields gene and pLI 
  extracted from the file at 
    
    ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/functional_gene_constraint/
      fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt

  To use another values file, add it as a parameter i.e.

     ./vep -i variants.vcf --plugin pLI,values_file.txt
    


=cut

package pLI;

use strict;
use warnings;

use DBI;
use Data::Dumper;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use List::MoreUtils qw/zip/;


my %include_columns = (
  "gene" => {
    "name" => "pLI_gene value"
  }
);
sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);
  
  my $file = $self->params->[0];
  
  my %scores;

  
  my $value = $self->params->[1] if (defined ($self->params->[1]));
  
  if(!$file) {
    my $plugin_dir = $INC{'pLI.pm'};
    $plugin_dir =~ s/pLI\.pm//i;
    $file = $plugin_dir.'/pLI_values.txt';
  }
  
  die("ERROR: pLI values file $file not found\n") unless $file && -e $file;

  # to get only the first line
  open IN, "<",  $file;
  while (<IN>){
    next unless  m/gene/;
    chomp;
    $_ =~   m/gene|pLI/;
    $self->{headers} = [split];
    
  }
 
  close IN;


  die "ERROR: Could not read headers from $file\n" unless defined($self->{headers});


  die "Error: File does not have a gene column " unless grep {$_ eq "gene"} @{$self->{headers}};
  $self->{header}{$include_columns{"gene"}{"name"}}  =  "pLI value by gene";
  open my $fh, "<", $file;
  while (<$fh>) {
    chomp;
    my ($gene, $score) = split;
    next if $score eq 'pLI';
    $scores{lc($gene)} = sprintf("%.2f", $score);
  } 
  close $fh;
  
  
  die("ERROR: No scores read from $file\n") unless scalar keys %scores;

  $self->{scores} = \%scores;
 
 
  return $self;


  
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  return $self->{header};
}

sub run {
  my $self = shift;
  my $tva = shift;
  

  my $symbol = $tva->transcript->{_gene_symbol} || $tva->transcript->{_gene_hgnc};
  return {} unless $symbol;
  return $self->{scores}->{lc($symbol)} ? { $include_columns{"gene"}{"name"} => $self->{scores}->{lc($symbol)}} : {};
  
}



1;

