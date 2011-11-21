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
    return ['Motif', 'Transcript', 'RegulatoryFeature'];
}

sub run {
    my ($self, $mfva) = @_;

    return undef unless $mfva->isa('Bio::EnsEMBL::Variation::MotifFeatureVariationAllele');

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

