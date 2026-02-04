=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2026] EMBL-European Bioinformatics Institute

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

 NearestExonJB

=head1 SYNOPSIS

 mv NearestExonJB.pm ~/.vep/Plugins
 ./vep -i variations.vcf --cache --plugin NearestExonJB

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 finds the nearest exon junction boundary to a coding sequence variant. More than 
 one boundary may be reported if the boundaries are equidistant or if using option
 --intronic.

 The plugin will report the Ensembl identifier of the exon, the distance to the
 exon boundary, the boundary type (start or end of exon) and the total
 length in nucleotides of the exon.

 Various key=value parameters can be altered by passing them to the plugin command:

   max_range : maximum search range in bp (default: 10000)

   intronic  : set to 1 to check nearest exons for intronic variants (default: 0)
               returns the nearest exon upstream and downstream without considering
               the max_range.


 Parameters are passed e.g.:

   --plugin NearestExonJB,max_range=50000
   --plugin NearestExonJB,max_range=50000,intronic=1
   --plugin NearestExonJB,intronic=1

=cut

package NearestExonJB;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my $char_sep = "|";

my %CONFIG = (
  max_range => 10000,
  intronic => 0
);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $params = $self->params;

  # get output format
  $char_sep = "+" if ($self->{config}->{output_format} eq 'vcf');

  foreach my $param(@$params) {
    my ($key, $val) = split('=', $param);
    die("ERROR: Failed to parse parameter $param\n") unless defined($key) && defined($val);
    $CONFIG{$key} = $val;
  }

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $header;

  if($CONFIG{intronic} == 1) {
    $header = 'Nearest Exon Junction Boundary. Format:';
  }
  else {
    $header = 'Nearest Exon Junction Boundary (coding sequence variants only). Format:';
  }

  $header .= join($char_sep, qw(ExonID distance start/end length) );

  return {
    NearestExonJB => $header,
  }
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->base_variation_feature;
  my $trv = $tva->base_transcript_variation;

  my $loc_string = sprintf("%s:%s-%i-%i", $trv->transcript_stable_id, $vf->{chr} || $vf->seq_region_name, $vf->{start}, $vf->{end});

  if(!exists($self->{_cache}) || !exists($self->{_cache}->{$loc_string})) {
    my $exons = $trv->_overlapped_exons; # intronic variants do not overlap any exon
    my %dists;
    my $min = $CONFIG{max_range};

    # For option --intronic, fetch the list of exons with different method
    # Do not take into account the max_range
    if(scalar @{$exons} == 0 && $CONFIG{intronic} == 1) {
      my $intron_numbers = $trv->intron_number();
      my $consequences = join(",", map { $_->SO_term } @{$tva->get_all_OverlapConsequences});

      if(defined $intron_numbers && $consequences =~ /intron/) {
        $exons = $trv->_sorted_exons;
        my ($intron_number, $total_number) = split(/\//, $intron_numbers);
 
        # Get the number of exons before and after the intron
        my $exon_before = $intron_number;
        my $exon_after = $intron_number + 1;

        my @exons_tmp;
        # In the reverse strand we get the last two exons from the list
        if($tva->transcript->strand < 0) {
          push(@exons_tmp, $exons->[-$exon_before]);
          push(@exons_tmp, $exons->[-$exon_after]);
        }
        else {
          push(@exons_tmp, $exons->[$exon_before -1]);
          push(@exons_tmp, $exons->[$exon_after -1]);
        }

        $exons = \@exons_tmp;
      }
    }

    foreach my $exon (@$exons) {
      my $startD = abs ($vf->start - $exon->seq_region_start);
      my $endD = abs ($vf->end - $exon->seq_region_end);
      if ($startD < $endD){
        $dists{$exon->stable_id}{$startD} = 'start';
        $dists{$exon->stable_id}{len} = $exon->length;
        $min = $startD if $min > $startD;
      } elsif ($startD > $endD){
        $dists{$exon->stable_id}{$endD} = 'end';
        $dists{$exon->stable_id}{len} = $exon->length;
        $min = $endD if $min > $endD;
      } else {
        $dists{$exon->stable_id}{$startD} = "start_end";
        $dists{$exon->stable_id}{len} = $exon->length;
        $min = $startD if $min > $startD;
      }
    }

    my @finalRes;
    # For option --intronic, return the closest exons (upstream/dowsntream) from the intron
    if(scalar @{$exons} == 2 && $CONFIG{intronic} == 1) {
      foreach my $exon (keys %dists) {
        my $inner_hash = $dists{$exon};
        my $length_value;
        my $type;
        my $distance_value;

        for my $internal_key (keys %{$inner_hash}) {
          if($internal_key eq "len") {
            $length_value = $inner_hash->{$internal_key};
          }
          else {
            $type = $inner_hash->{$internal_key};
            $distance_value = $internal_key;
          }
        }

        my $string = $exon . $char_sep . $distance_value . $char_sep . $type . $char_sep . $length_value;
        push(@finalRes, $string);
      }
    }
    else {
        # This is the default behaviour of the plugin
        foreach my $exon (keys %dists){
        if (exists $dists{$exon}{$min}) {
          push(@finalRes, $exon.$char_sep.$min.$char_sep.$dists{$exon}{$min}.$char_sep.$dists{$exon}{len})
        }
      }
    }

    $self->{_cache}->{$loc_string} = scalar @finalRes ? join(",", @finalRes) : undef;
  }
  return $self->{_cache}->{$loc_string} ? { NearestExonJB => $self->{_cache}->{$loc_string} } : {};
}

1;
