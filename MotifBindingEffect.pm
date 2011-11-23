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

 MotifBindingEffect

=head1 SYNOPSIS

 mv MotifBindingEffect.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin MotifBindingEffect

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that calculates
 the difference in binding affinity for the reference and variant sequences for
 a variant that overlaps a MotifFeature, Ensembl's representation of a TF binding
 motif mapped to a particular genomic location. It adds two new entry classes to 
 the VEP's Extra column, BINDING_SCORE_DELTA which is the change in relative 
 binding affinity, and MOTIF_POSITION which gives the relative position of the 
 variant in the motif.

=cut

package MotifBindingEffect;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub version {
    return '2.3';
}

sub get_header_info {
    return {
        BINDING_SCORE_DELTA => "The difference in motif binding affinity between the reference and alternative sequences",
        MOTIF_POSITION      => "The relative position of the variation in the motif"
    };
}

sub feature_types {
    return ['Motif'];
}

sub run {
    my ($self, $mfva) = @_;

    my $vf = $mfva->motif_feature_variation->variation_feature;
    my $mf = $mfva->motif_feature;

    my $allele_seq      = $mfva->feature_seq;
    my $ref_allele_seq  = $mfva->motif_feature_variation->get_reference_MotifFeatureVariationAllele->feature_seq;

    if ($allele_seq eq '-' || 
        $ref_allele_seq eq '-' || 
        length($allele_seq) != length($ref_allele_seq)) {
        # we can't call a score because the sequence will change length
        return undef;
    }

    # get the 0-based relative start of the vf with respect to this mf

    my $mf_start = $vf->seq_region_start - $mf->seq_region_start;

    # adjust if the motif is on the reverse strand

    $mf_start = $mf->binding_matrix->length - $mf_start - 1 if $mf->strand < 0;

    return {} if $mf_start < 0;

    my $var_len = length($allele_seq);

    return {} if $var_len > $mf->length;

    my $mf_seq = $mf->seq;

    my $matrix = $mf->binding_matrix;

    # get the binding affinity of the reference sequence
    my $ref_affinity = $matrix->relative_affinity($mf_seq);

    # splice in the variant sequence
    substr($mf_seq, $mf_start, $var_len) = $allele_seq;

    # and get the affinity of the variant sequence
    my $var_affinity = $matrix->relative_affinity($mf_seq);

    # we report the motif position in 1-based coordinates

    my $delta = ($var_affinity - $ref_affinity);

    return {
        MOTIF_POSITION      => $mf_start + 1,
        BINDING_SCORE_DELTA => sprintf "%.3f", $delta
    };
}

1;

