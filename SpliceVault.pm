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

 SpliceVault

=head1 SYNOPSIS

 mv SpliceVault.pm ~/.vep/Plugins

 ./vep -i variations.vcf --plugin SpliceVault,file=/path/to/SpliceVault_data.tsv.gz

 # Stringely select predicted loss-of-function (pLoF) splicing variants
 ./filter_vep -i variant_effect_output.txt --filter "SPLICEVAULT_OUT_OF_FRAME_EVENTS >= 3"

=head1 DESCRIPTION

 A VEP plugin that retrieves SpliceVault data to predict exon-skipping events
 and activated cryptic splice sites based on the most common mis-splicing events
 around a splice site.

 This plugin returns the most common (top 4) variant-associated mis-splicing
 events for variants that overlap a near-splice-site region (-17:+3 for splice
 acceptor sites and -3:+8 for splice donor sites) based on SpliceVault data.
 Each event includes the following information:
   - Type: exon skipping (ES), cryptic donor (CD) or cryptic acceptor (CA)
   - Description: exons skiped (ES events) or cryptic distance to the annotated
     splice site (CD/CA events)
   - Percent of supporting samples: samples supporting the event over total
     samples where splicing occurs in that site (in percentage)
   - Frameshift: inframe or out-of-frame event

 Please cite the SpliceVault publication alongside the VEP if you use this
 resource: https://pubmed.ncbi.nlm.nih.gov/36747048/

 The tabix utility must be installed in your path to use this plugin. The
 SpliceVault TSV and respective index (TBI) for GRCh38 can be downloaded from
 https://ftp.ensembl.org/pub/current_variation/SpliceVault/SpliceVault_data.tsv.gz

 To filter results, please use filter_vep with the output file or standard
 output. Documentation on filter_vep is available at:
 https://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html

=cut

package SpliceVault;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();
  
  # Check files in arguments
  my @files; 
  my $params = $self->params_to_hash();

  for my $key (keys %{$params}) {
    push @files, $params->{$key};
  }

  die "ERROR: add Tabix-indexed file, for example:\n" .
    "  --plugin SpliceVault,file=/path/to/SpliceVault_data.tsv.gz \n"
    unless @files > 0;
  $self->add_file($_) for @files;

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {
    'SpliceVault_sample_count' => 'Number of SpliceVault annotated splicing samples',
    'SpliceVault_out_of_frame_events' => 'Number of top SpliceVault events that are out of frame',
    'SpliceVault_top*_event' => 'SpliceVault top events with the following colon-separated fields: type:description:percent_of_supporting_samples:event_frame',
  }
}

sub _join_results {
  my $self = shift;
  my $all_results = shift;
  my $res = shift;

  if ($self->{config}->{output_format} eq 'json' || $self->{config}->{rest}) {
    for (keys %$res) {
      next unless $_ =~ /TOP[0-9]+_EVENT/;
      my @fields = split /:/, $res->{$_};
      delete $res->{$_};
      $res->{$_}->{'TYPE'} = $fields[0];
      $res->{$_}->{'DESCRIPTION'} = $fields[1];
      $res->{$_}->{'PERCENT_OF_SUPPORTING_SAMPLES'} = $fields[2];
      $res->{$_}->{'FRAMESHIFT'} = $fields[3];
    }
    # Group results for JSON and REST
    $all_results->{"SpliceVault"} = [] unless defined $all_results->{"SpliceVault"};
    push(@{ $all_results->{"SpliceVault"} }, $res);
  } else {
    # Create array of results per key
    for (keys %$res) {
      $all_results->{$_} = [] if !$all_results->{$_};
      push(@{ $all_results->{$_} }, $res->{$_} || "NA");
    }
  }
  return $all_results;
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  my @data = @{ $self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end}) };

  my $all = {};
  foreach (@data) {
    next unless $tva->transcript->stable_id eq $_->{feature};
    $all = $self->_join_results($all, $_->{result});
  }
  return $all;
}

sub parse_data {
  my ($self, $line) = @_;
  my ($feature, $chrom, $start, $end, $count, $out_of_frame, @top_events) = split /\t/, $line;

  my $res = {
    feature => $feature,
    result  => {
      SpliceVault_sample_count        => $count,
      SpliceVault_out_of_frame_events => $out_of_frame,
    }
  };

  my $n = 0;
  for my $event (@top_events) {
    $n++;
    $event =~ s/;/:/g;
    
    if ($event =~ /^ES/) {
      $event =~ s/^ES/exon_skipping/g;
    } elsif ($event =~ /^CA/) {
      $event =~ s/^CA/cryptic_acceptor/g;
    } elsif ($event =~ /^CD/) {
      $event =~ s/^CD/cryptic_donor/g;
    }
    $res->{result}->{"SpliceVault_top${n}_event"} = $event;
  }
  return $res;
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
