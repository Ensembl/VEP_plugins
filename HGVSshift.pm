=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

  HGVSshift

=head1 SYNOPSIS

  mv HGVSshift.pm ~/.vep/Plugins
  ./vep -i variations.vcf --cache --plugin HGVSshift

=head1 DESCRIPTION

  This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
  reports the "opposite" 3'-shifted HGVS notation to the main HGVSc and
  HGVSp VEP notations, i.e.

  a) If using "--shift_hgvs 1" (default), "unshifted" HGVS notations
     will be reported as HGVSc_unshifted and HGVSp_unshifted.

  b) If using "--shift_hgvs 0", "shifted" HGVS notations will be
     reported as HGVSc_shifted and HGVSp_shifted

  New notations are only reported if they differ from the main notations, and
  are reported in addition to, not instead of, the main notations.

=cut

package HGVSshift;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::TranscriptVariationAllele;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  
  # check config is OK
  
  # FASTA file defined, optimal
  if(!defined($self->{config}->{fasta})) {
    
    # offline mode won't work without FASTA
    die("ERROR: Cannot generate HGVS without either a FASTA file (--fasta) or a database connection (--cache or --database)\n") if defined($self->{config}->{offline}) and !defined($self->{config}->{quiet});
    
    # cache mode will work, but DB will be accessed
    warn("WARNING: Database will be accessed using this plugin; use a FASTA file (--fasta) for optimal performance\n") if defined($self->{config}->{cache}) and !defined($self->{config}->{quiet});
  }
  
  #
  if(!defined($self->{config}->{hgvs})) {
    warn("WARNING: Plugin is enabling --hgvs\n") unless defined($self->{config}->{quiet});
    $self->{config}->{hgvs} = 1;
  }
  
  $self->{remove_transcript_ID} = $self->params->[0];
  
  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub variant_feature_types {
  return ['VariationFeature'];
}

sub get_header_info {
  my $self = shift;
  my $prefix = $self->shifting_enabled ? 'Uns' : 'S';
  
  return {
    'HGVSc_'.lc($prefix).'hifted' => $prefix.'hifted HGVS transcript notation',
    'HGVSp_'.lc($prefix).'hifted' => $prefix.'hifted HGVS protein notation',
  };
}

sub run {
  my ($self, $tva) = @_;
  
  # check var class, shifting only happens to insertions and deletions
  my $var_class = $tva->variation_feature->var_class();
  return {} unless $var_class eq 'insertion' || $var_class eq 'deletion';
  
  # find out config status
  my $shifting_enabled = $self->shifting_enabled;
  
  my $hgvs = {};
  
  # shifting enabled (default from e!80 onwards)
  if($shifting_enabled) {
    
    # we only want to report unshifted if this one has been shifted
    if($tva->can('hgvs_offset') && $tva->hgvs_offset) {
      $self->reset_hgvs($tva);
      $self->switch_shifting_state($tva, 0);
      $hgvs->{HGVSc_unshifted} = $tva->hgvs_transcript;
      $hgvs->{HGVSp_unshifted} = $tva->hgvs_protein;
      $self->switch_shifting_state($tva, 1);
    }
  }
  
  # shifting disabled
  else {
    
    # get the original ones so we don't report the same thing twice
    # unfortunately we do have to calculate it twice...
    my ($original_hgvsc, $original_hgvsp) = ($tva->hgvs_transcript, $tva->hgvs_protein);
    
    $self->reset_hgvs($tva);
    $self->switch_shifting_state($tva, 1);
    
    my ($new_hgvsc, $new_hgvsp) = ($tva->hgvs_transcript, $tva->hgvs_protein);
    $hgvs->{HGVSc_shifted} = $new_hgvsc if $new_hgvsc && $original_hgvsc && $new_hgvsc ne $original_hgvsc;
    $hgvs->{HGVSp_shifted} = $new_hgvsp if $new_hgvsp && $original_hgvsp && $new_hgvsp ne $original_hgvsp;
    $self->switch_shifting_state($tva, 0);    
  }
  
  # delete empty keys
  delete $hgvs->{$_} for grep {!$hgvs->{$_}} keys %$hgvs;
  
  return $hgvs;
}

sub shifting_enabled {
  my $self = shift;
  
  if(!exists($self->{shifting_enabled})) {
    my $shifting_enabled;
    my $config = $self->{config};
  
    if(defined($config->{shift_hgvs})) {
      $shifting_enabled = $config->{shift_hgvs};
    }
    elsif(defined($Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME)) {
      $shifting_enabled = $Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME;
    }
    elsif(defined($Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME)) {
      $shifting_enabled = $Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME;
    }
    
    $self->{shifting_enabled} = $shifting_enabled;
  }
  
  return $self->{shifting_enabled};
}

sub switch_shifting_state {
  my $self = shift;
  my $tva = shift;
  my $newval = shift;
  
  my $tv = $tva->transcript_variation;
  
  $tv->adaptor->db->shift_hgvs_variants_3prime($newval) if defined $tv->adaptor() && UNIVERSAL::can($tv->adaptor, 'isa');
  
  no warnings 'once';
  $Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME = $newval;
  no warnings 'once';
  $Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME = $newval;
  
  
}

sub reset_hgvs {
  my $self = shift;
  my $tva = shift;
  delete $tva->{$_} for grep {/^hgvs/} keys %$tva;
}

1;
