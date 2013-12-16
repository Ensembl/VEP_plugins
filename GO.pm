=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and                                                   
 Genome Research Limited.  All rights reserved.                                                                      

 This software is distributed under a modified Apache license.                                                       
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html                                                               

=head1 CONTACT                                                                                                       

 Will McLaren <wm2@ebi.ac.uk>
    
=cut

=head1 NAME

 GO

=head1 SYNOPSIS

 mv GO.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin GO

=head1 DESCRIPTION

 A VEP plugin that retrieves Gene Ontology terms associated with
 transcripts/translations via the Ensembl API. Requires database connection.
 
=cut

package GO;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  
  # connect to DB for offline users
  my $config = $self->{config};
  my $reg = $config->{reg};
  
  if(!defined($self->{config}->{sa})) {
    $reg = 'Bio::EnsEMBL::Registry';
    $reg->load_registry_from_db(
      -host       => $config->{host},
      -user       => $config->{user},
      -pass       => $config->{password},
      -port       => $config->{port},
      -db_version => $config->{db_version},
      -species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
      -verbose    => $config->{verbose},
      -no_cache   => $config->{no_slice_cache},
    );
  }
  
  return $self;
}

sub version {
  return 73;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return { 'GO' => 'GO terms associated with protein product'};
}

sub run {
  my ($self, $tva) = @_;
  
  my $tr = $tva->transcript->translation;
  return {} unless defined($tr);
  
  my $entries = $tr->get_all_DBEntries('GO');
  
  my $string = join(",", map {$_->display_id.':'.$_->description} @$entries);
  $string =~ s/\s+/\_/g;
  
  return { GO => $string };
}

1;

