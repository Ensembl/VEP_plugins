=head1 LICENSE
                                                                                                                     
 Copyright (c) 1999-2011 The European Bioinformatics Institute and                                                   
 Genome Research Limited.  All rights reserved.                                                                      
                                                                                                                     
 This software is distributed under a modified Apache license.                                                       
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html                                                               
                                                                                                                     
=head1 CONTACT                                                                                                       

 Graham Ritchie <grsr@ebi.ac.uk>
    
=cut

=head1 NAME

 CCDSFilter

=head1 SYNOPSIS

 mv CCDSFilter.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin CCDSFilter

=head1 DESCRIPTION

 A simple VEP filter plugin that limits output to variants that
 fall in transcripts which have CCDS coding sequences.

=cut

package CCDSFilter;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepFilterPlugin);

sub feature_types {
    return ['Transcript'];
}

sub include_line {
    my ($self, $tva) = @_;

    my $t = $tva->transcript;

    my @entries = grep {$_->database eq 'CCDS'} @{$t->get_all_DBEntries};
   
    return scalar @entries;
}

1;

