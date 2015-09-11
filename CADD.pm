=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

 Will McLaren <wm2@ebi.ac.uk>
    
=cut

=head1 NAME

 CADD

=head1 SYNOPSIS

 mv CADD.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz

=head1 DESCRIPTION

 A VEP plugin that retrieves CADD scores for variants from one or more
 tabix-indexed CADD data files.
 
 Please cite the CADD publication alongside the VEP if you use this resource:
 http://www.ncbi.nlm.nih.gov/pubmed/24487276
 
 The tabix utility must be installed in your path to use this plugin. The CADD
 data files can be downloaded from
 http://cadd.gs.washington.edu/download
 
=cut

package CADD;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  
  # test tabix
  die "ERROR: tabix does not seem to be in your path\n" unless `which tabix 2>&1` =~ /tabix$/;
  
  # get CADD file
  my $files = $self->params;
  
  foreach my $file(@$files) {
    
    # remote files?
    if($file =~ /tp\:\/\//) {
      my $remote_test = `tabix -f $file 1:1-1 2>&1`;
      if($remote_test && $remote_test !~ /get_local_version/) {
        die "$remote_test\nERROR: Could not find file or index file for remote annotation file $file\n";
      }
    }

    # check files exist
    else {
      die "ERROR: CADD file $file not found\n" unless -e $file;
      die "ERROR: Tabix index file $file\.tbi not found - perhaps you need to create it first?\n" unless -e $file.'.tbi';
    }
  }
  
  $self->{files} = $files;
  
  return $self;
}

sub feature_types {
  return ['Feature','Intergenic'];
}

sub get_header_info {
  my $self = shift;
  return {
    CADD_PHRED => 'PHRED-like scaled CADD score',
    CADD_RAW   => 'Raw CADD score'
  }
}

sub run {
  my ($self, $tva) = @_;
  
  my $vf = $tva->variation_feature;
  
  # get allele, reverse comp if needed
  my $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;
  
  return {} unless $allele =~ /^[ACGT-]$/;
  
  # adjust coords to account for VCF-like storage of indels
  my ($s, $e) = ($vf->{start} - 1, $vf->{end} + 1);
  
  my $pos_string = sprintf("%s:%i-%i", $vf->{chr}, $s, $e);
  
  # clear cache if it looks like the coords are the same
  # but allele type is different
  delete $self->{cache} if
    defined($self->{cache}->{$pos_string}) &&
    scalar keys %{$self->{cache}->{$pos_string}} &&
    !defined($self->{cache}->{$pos_string}->{$allele});
  
  my %cadd_data;
  
  # cached?
  if(defined($self->{cache}) && defined($self->{cache}->{$pos_string})) {
    %cadd_data = %{$self->{cache}->{$pos_string}};
  }
  
  # read from file(s)
  else {
    foreach my $file(@{$self->{files}}) {
      open TABIX, sprintf("tabix -f %s %s |", $file, $pos_string);
    
      while(<TABIX>) {
        chomp;
        s/\r$//g;
        my ($c, $s, $ref, $alt, $raw, $phred) = split /\t/;
      
        # do VCF-like coord adjustment for mismatched subs
        my $e = ($s + length($ref)) - 1;
        if(length($alt) != length($ref)) {
          $s++;
          $ref = substr($ref, 1);
          $alt = substr($alt, 1);
          $ref ||= '-';
          $alt ||= '-';
        }
      
        next unless $s == $vf->{start} && $e == $vf->{end};
      
        $cadd_data{$alt} = {
          CADD_RAW   => $raw,
          CADD_PHRED => $phred
        };
      }
    
      close TABIX;
      
      last if scalar keys %cadd_data;
    }
  }
  
  # overwrite cache
  $self->{cache} = {$pos_string => \%cadd_data};
  
  return defined($cadd_data{$allele}) ? $cadd_data{$allele} : {};
}

1;

