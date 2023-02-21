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

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>
    
=cut

=head1 NAME

 TranscriptAnnotator.pm

=head1 SYNOPSIS

 mv TranscriptAnnotator.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin TranscriptAnnotator,/path/to/file.txt.gz

=head1 DESCRIPTION

 A VEP plugin that annotates transcript consequences based on a given file:
     
   --plugin TranscriptAnnotator,${HOME}/file.txt.gz
 
 The tabix and bgzip utilities must be installed in your path to read the
 annotation: check https://github.com/samtools/htslib.git for installation
 instructions.
 
=cut

package TranscriptAnnotator;

use strict;
use warnings;
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);
use Data::Dumper;

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);
  my $file = $self->params->[0]; 
  $self->add_file($file);
  $self->get_user_params();
  return $self;
}

sub version {
  return 106;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {
      SIFT_score => 'SIFT score',
      SIFT_prediction => 'SIFT prediction'
  };
}

sub run {
  my ($self, $tva) = @_;
  my $tr            = $tva->transcript;
  my $transcript_id = $tr->{stable_id};

  my $vf            = $tva->variation_feature;
  my $seqname       = $vf->{slice}->{seq_region_name};
  my $ref           = $vf->ref_allele_string;
  my ($alt,   $len) = map {$_ => 1} @{$vf->alt_alleles};
  my ($start, $end) = ($vf->{start}, $vf->{end});
  
  my @data = @{$self->get_data($seqname, $start, $end)};
  foreach (@data) {
    return $_->{result} if $_->{seqname} eq $seqname &&
                           $_->{start} eq $start &&
                           $_->{end} eq $end &&
                           $_->{ref} eq $ref &&
                           $_->{alt} eq $alt &&
                           $_->{transcript_id} eq $transcript_id;
  }
  return {};
}

sub _trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };

sub parse_data {
  my ($self, $line) = @_;
  my ($seqname, $pos, $ref, $alt, $transcript_id, $source, $score, $pred) = split /\t/, $line;

  return {
    seqname => _trim($seqname),
    start => _trim($pos),
    end => _trim($pos),
    ref => _trim($ref),
    alt => _trim($alt),
    transcript_id => _trim($transcript_id),
    result => {
      SIFT_score => _trim($score),
      SIFT_prediction => _trim($pred)
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
