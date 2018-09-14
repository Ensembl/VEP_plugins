=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

 Stephen Kazakoff <sh.kazakoff@gmail.com>
    
=cut

=head1 NAME

 LRG_RefSeqGene

=head1 SYNOPSIS

 mv LRG.pm ~/.vep/Plugins
 ./vep -i variants.vcf --plugin LRG_RefSeqGene,/path/to/LRG_RefSeqGene

=head1 DESCRIPTION

 A VEP plugin that retrieves the LRG_RefSeqGene category from the
 LRG_RefSeqGene file. You can obtain the LRG_RefSeqGene file using:

 > wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene

 For more information about RefSeqGene and LRG (Locus Reference Genomic),
 please see:

 https://www.ncbi.nlm.nih.gov/refseq/rsg/lrg

=cut

package LRG_RefSeqGene;

use strict;
use warnings;

use Text::CSV;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $params = $self->params;

  my $file = shift @$params;

  die("ERROR: LRG_RefSeqGene $file not found\n") unless $file && -e $file;

  $self->{file} = $file;

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {
    LRG_RefSeqGene => "LRG_RefSeqGene category",
  };
}

sub run {
  my ($self, $tva) = @_;

  my $tr = $tva->transcript;

  my $gene_id = $tr->{_gene_stable_id} || $tr->{_gene}->stable_id;
  my $transcript_id = $tr->stable_id;
  my $protein_id = $tr->{_protein} || '';

  unless (defined($self->{cache})) {

    open(my $fh, $self->{file}) or die $!;
    my $tsv = Text::CSV->new({ binary => 1, auto_diag => 1, sep_char => "\t" });
    $tsv->column_names($tsv->getline($fh));

    my %data;

    while (my $href = $tsv->getline_hr($fh)) {
      my ($gene_id, $transcript_id, $protein_id, $category) = @{$href}{qw(GeneID RNA Protein Category)};

      $data{$gene_id,$transcript_id,$protein_id}{$category} = undef;
    }
    close $fh;

    $self->{cache} = \%data;
  }

  if (defined(my $categories = $self->{cache}{$gene_id,$transcript_id,$protein_id})) {

    return { LRG_RefSeqGene => join(',', sort keys %$categories) };
  }

  return {};
}

1;

