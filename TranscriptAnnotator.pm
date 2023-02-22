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

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);
  my $file = $self->params->[0]; 
  $self->add_file($file);
  $self->get_user_params();

  return $self;
}

sub version {
  return 110;
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
  my $tr     = $tva->transcript;

  # Get transcript ID for Ensembl and RefSeq
  my @refseq = split(/,/, $tr->{_refseq}) unless $tr->{_refseq} eq '-';
  my @tr_id  = ( $tr->{stable_id}, @refseq );
  my $vf     = $tva->variation_feature;

  # Get allele
  my $alt_allele = $tva->variation_feature_seq;
  my $ref_allele = $vf->ref_allele_string;

  my @data = @{ $self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end}) };
  return {} unless(@data);

  foreach my $var (@data) {
    my $is_same_transcript = grep { $var->{transcript_id} eq $_ } @tr_id;

    my $matches = get_matched_variant_alleles(
       {
         ref    => $ref_allele,
         alts   => [$alt_allele],
         pos    => $vf->{start},
         strand => $vf->strand
       },
       {
        ref  => $var->{ref},
        alts => [$var->{alt}],
        pos  => $var->{start},
       }
     );

    return $var->{result} if $is_same_transcript && (@$matches);
  }
  return {};
}

sub _trim_whitespaces { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };

sub parse_data {
  my ($self, $line) = @_;

  my @data = split /\t/, $line;
  @data = map(_trim_whitespaces($_), @data);

  my ($seqname, $pos, $ref, $alt, $transcript_id, $source, $score, $pred) = @data;

  return {
    seqname => $seqname,
    start => $pos,
    end => $pos,
    ref => $ref,
    alt => $alt,
    transcript_id => $transcript_id,
    result => {
      SIFT_score => $score,
      SIFT_prediction => $pred
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
