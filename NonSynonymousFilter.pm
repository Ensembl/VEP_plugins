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

 NonSynonymousFilter

=head1 SYNOPSIS

 mv NonSynonymousFilter.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin NonSynonymousFilter

=head1 DESCRIPTION

 A simple example VEP filter plugin that limits output to non-synonymous variants

=cut

package NonSynonymousFilter;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepFilterPlugin);

sub feature_types {
    return ['Transcript'];
}

sub include_line {
    my ($self, $vfoa) = @_;

    return 0 unless $vfoa->can('pep_allele_string');

    if (my $pep_alleles = $vfoa->pep_allele_string) {
        return $pep_alleles =~ /\//;
    }

    return 0;
}

1;

