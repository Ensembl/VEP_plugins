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

 TSSDistance

=head1 SYNOPSIS

 mv TSSDistance.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin CCDSFilter

=head1 DESCRIPTION

 A VEP plugin that calculates the distance from the transcription
 start site for upstream variants.

=cut

package TSSDistance;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub get_header_info {
    return {
        TSSDistance => "Distance from the transcription start site"
    };
}

sub feature_types {
    return ['Transcript'];
}

sub run {
    my ($self, $tva) = @_;

    my $t = $tva->transcript;
    my $vf = $tva->variation_feature;

    my $dist;

    if ($t->strand == 1) {
        $dist = $t->start - $vf->end;
    }
    else {
        $dist = $vf->start - $t->end;
    }

    if ($dist > 0) {
        return {
            TSSDistance => $dist,
        }
    }
    else {
        return {};
    }
}

1;
