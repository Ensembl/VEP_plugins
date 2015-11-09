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

 FATHMM_MKL

=head1 SYNOPSIS

 mv FATHMM_MKL.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i input.vcf --plugin FATHMM_MKL,fathmm-MKL_Current.tab.gz

=head1 DESCRIPTION

 A VEP plugin that retrieves FATHMM-MKL scores for variants from a tabix-indexed
 FATHMM-MKL data file.
 
 See https://github.com/HAShihab/fathmm-MKL for details.

 NB: The currently available data file is for GRCh37 only.
 
=cut

package FATHMM_MKL;

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
  
  # get FATHMM_MKL file
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
      die "ERROR: FATHMM-MKL file $file not found\n" unless -e $file;
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
    FATHMM_MKL_C => 'FATHMM-MKL coding score',
    FATHMM_MKL_NC => 'FATHMM-MKL non-coding score',
  }
}

sub run {
  my ($self, $tva) = @_;
  
  my $vf = $tva->variation_feature;

  return {} unless $vf->{start} eq $vf->{end};
  
  # get allele, reverse comp if needed
  my $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;
  
  return {} unless $allele =~ /^[ACGT]$/;
  
  # adjust coords, file is BED-like (but not 0-indexed, go figure...)
  my ($s, $e) = ($vf->{start}, $vf->{end} + 1);
  
  my $pos_string = sprintf("%s:%i-%i", $vf->{chr}, $s, $e);
  
  # clear cache if it looks like the coords are the same
  # but allele type is different
  delete $self->{cache} if
    defined($self->{cache}->{$pos_string}) &&
    scalar keys %{$self->{cache}->{$pos_string}} &&
    !defined($self->{cache}->{$pos_string}->{$allele});
  
  my %fathmm_data;
  
  # cached?
  if(defined($self->{cache}) && defined($self->{cache}->{$pos_string})) {
    %fathmm_data = %{$self->{cache}->{$pos_string}};
  }
  
  # read from file(s)
  else {
    foreach my $file(@{$self->{files}}) {
      open TABIX, sprintf("tabix -f %s %s |", $file, $pos_string);
    
      while(<TABIX>) {
        chomp;
        s/\r$//g;
        my ($c, $s, $e, $ref, $alt, $nc_score, $nc_groups, $c_score, $c_groups) = split /\t/;

        # BED-like adjustment
        $e--;
      
        next unless $s == $vf->{start} && $e == $vf->{end};
      
        $fathmm_data{$alt} = {
          FATHMM_MKL_C  => $c_score,
          FATHMM_MKL_NC => $nc_score,
        };
      }
    
      close TABIX;
      
      last if scalar keys %fathmm_data;
    }
  }
  
  # overwrite cache
  $self->{cache} = {$pos_string => \%fathmm_data};
  
  return defined($fathmm_data{$allele}) ? $fathmm_data{$allele} : {};
}

1;

