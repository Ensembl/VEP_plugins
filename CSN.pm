=head1 LICENSE
                                                                                                                     
 Copyright (c) 1999-2015 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.                                                                      
                                                                                                                     
 This software is distributed under a modified Apache license.                                                       
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html                                                               
                                                                                                                     
=head1 CONTACT                                                                                                       

 Will McLaren <wm2@ebi.ac.uk>
    
=cut

=head1 NAME

  CSN

=head1 SYNOPSIS

  mv CSN.pm ~/.vep/Plugins
  perl variant_effect_predictor.pl -i variations.vcf --cache --plugin CSN

=head1 DESCRIPTION

  This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
  reports Clinical Sequencing Nomenclature (CSN) for variants.

  Each notation is given with reference to the transcript identifier;
  specify "--plugin CSN,1" to remove this identifier from the CSN string.

  You may also wish to specify "--no_escape" to prevent the "=" in "p.="
  notations being converted to the URI-escaped equivalent "p.%3D"; doing
  so may break parsers looking for "=" as a KEY=VALUE separator.

  See http://biorxiv.org/content/early/2015/03/21/016808.1

=cut

package CSN;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  
  # check config is OK
  
  # FASTA file defined, optimal
  if(!defined($self->{config}->{fasta})) {
    
    # offline mode won't work without FASTA
    die("ERROR: Cannot generate CSN without either a FASTA file (--fasta) or a database connection (--cache or --database)\n") if defined($self->{config}->{offline});
    
    # cache mode will work, but DB will be accessed
    warn("WARNING: Database will be accessed using this plugin; use a FASTA file (--fasta) for optimal performance") if defined($self->{config}->{cache});
  }
  
  no warnings 'once';
  $Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME = 1;
  no warnings 'once';
  $Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME = 1;
  
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
  return { CSN => 'Clinical Sequencing Nomenclature'};
}

sub run {
  my ($self, $tva) = @_;
  
  my ($hgvs_c, $hgvs_p) = ($tva->hgvs_transcript || '', $tva->hgvs_protein || '');
  
  return {} unless $hgvs_c;
  
  # trim off transcript/protein ID
  $hgvs_c =~ s/.+\:// if $self->{remove_transcript_ID}; 
  $hgvs_p =~ s/.+\://;
  
  # change Ter to X
  $hgvs_p =~ s/Ter/X/g;
  
  # leave just p.=
  $hgvs_p = 'p.=' if $hgvs_p =~ /p\.\=/;
  
  # escape
  $hgvs_p =~ s/\=/\%3D/g unless $self->{config}->{no_escape};
  
  return { CSN => $hgvs_c.($hgvs_p ? '_'.$hgvs_p : '') };
}

1;
