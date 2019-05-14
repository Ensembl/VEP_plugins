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

 FunMotifs

=head1 SYNOPSIS

 mv FunMotifs.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin FunMotifs,/path/to/funmotifs/all_tissues.bed.gz,all,uterus
 ./vep -i variations.vcf --plugin FunMotifs,/path/to/funmotifs/blood.funmotifs_sorted.bed.gz,individual,fscore,dnase_seq
 
 Parameters Required:
 
 [0] : FunMotifs BED file
 [1] : Type of file provided ('all' or 'individual')
 [2]+ : List of columns to include within VEP output (e.g. fscore, skin, contactingdomain)
 

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds tissue-specific transcription factor motifs from FunMotifs to VEP output.

 Please cite the FunMotifs publication alongside the VEP if you use this resource, which you can find on:
 http://bioinf.icm.uu.se:3838/funmotifs/
 
 FunMotifs files can be downloaded from: http://bioinf.icm.uu.se:3838/funmotifs/

 The tabix utility must be installed in your path to use this plugin.

=cut
package FunMotifs;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);
  if(scalar(@{$self->params}) < 3){
    warn('Insufficient input parameters found for FunMotifs plugin');
    return $self;
  }
  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();

  my $file = $self->params->[0];

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return { FunMotifs => 'Annotated Transcription Factor Motifs'};
}

sub run {
  my ($self, $tva) = @_;

  my $params = $self->params;
  
  if(scalar(@{$self->params}) < 3){
    return $self;
  }
  
  my $filetype = ($params->[1] eq 'all') ? 'all' : 'individual';
  
  my $output_prefix = 'FM';
  
  $self->{motif_filetype} = $filetype;
  
  my $tissue_type = $params->[1];

  my $vf = $tva->variation_feature;

  my $column_headers = {
                         all => [qw(name score pval strand blood brain breast cervix colon esophagus kidney liver lung myeloid pancreas prostate skin stomach uterus)],
                         individual => [qw(name score pval strand fscore chromhmm contactingdomain dnase_seq fantom loopdomain numothertfbinding othertfbinding replidomain tfbinding tfexpr)],
                       };

  my ($res) = @{$self->get_data($vf->{chr}, $vf->{start}, $vf->{end})};

  shift $params; #Removes filename
  
  
  my %col_head_hash = map { $_ => 1 } @{$column_headers->{$filetype}};
  my $col_head_hashref = \%col_head_hash;
  my $output_hash = {};

  foreach my $para(@$params)
  {
      $output_hash->{$output_prefix . $para} = $res->{$para} if ($col_head_hashref->{$para} && defined($res->{$para}));
  } 

  return $output_hash;
}

sub parse_data {
  my ($self, $line) = @_;

  return $self->parse_data_all($line) if $self->{motif_filetype} eq 'all';
  return $self->parse_data_individual($line) if $self->{motif_filetype} eq 'individual';
}


sub parse_data_all {
  my ($self, $line) = @_;

  my ($chr, $start, $end, $name, $score, $strand, $pval, $blood, $brain, $breast, $cervix, $colon, $esophagus, $kidney, $liver, $lung, $myeloid, $pancreas, $prostate, $skin, $stomach, $uterus) = split /\t/, $line;

  return {
    chr => $chr,
    start => $start,
    end => $end,
    name => $name,
    score => $score,
    pval => $pval,
    strand => $strand,
    blood => $blood,
    brain => $brain,
    breast => $breast,
    cervix => $cervix,
    colon => $colon,
    esophagus => $esophagus,
    kidney => $kidney,
    liver => $liver,
    lung => $lung,
    myeloid => $myeloid,
    pancreas => $pancreas,
    prostate => $prostate,
    skin => $skin,
    stomach => $stomach,
    uterus  => $uterus,
  };
}


sub parse_data_individual {
  my ($self, $line) = @_;

  my ($chr, $start, $end, $name, $score, $strand, $pval, $fscore, $chromhmm, $contactingdomain, $dnase_seq, $fantom, $loopdomain, $numothertfbinding, $othertfbinding, $replidomain, $tfbinding, $tfexpr) = split /\t/, $line;

  return {
    name => $name,
    score => $score,
    pval => $pval,
    strand => $strand,
    fscore => $fscore,
    chromhmm => $chromhmm,
    contactingdomain => $contactingdomain,
    dnase_seq => $dnase_seq,
    fantom => $fantom,
    loopdomain => $loopdomain,
    numothertfbinding => $numothertfbinding,
    othertfbinding => $othertfbinding,
    replidomain => $replidomain,
    tfbinding => $tfbinding,
    tfexpr => $tfexpr,    
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
