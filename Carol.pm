package Carol;

use strict;
use warnings;

use Math::CDF qw(pnorm qnorm);

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my $CAROL_CUTOFF = 0.98;

sub version {
    return '2.3';
}

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {
        CAROL => "Combined Annotation scoRing toOL prediction",
    };
}

sub run {
    my ($self, $tva) = @_;
    
    my $pph_pred    = $tva->polyphen_prediction;
    my $pph_score   = $pph_pred ? ($pph_pred eq 'unknown' ? undef: $tva->polyphen_score) : undef;
    my $sift_score  = $tva->sift_score;

    my ($carol_pred, $carol_score) = compute_carol($pph_score, $sift_score);
    
    my $results = {};

    if (defined $carol_pred) {

        $carol_score = sprintf "%.3f", $carol_score;

        my $result = "$carol_pred ($carol_score)";

        if (@{ $self->params } > 0) {
            $result = $carol_pred if ($self->params->[0] =~ /^p/i);
            $result = $carol_score if ($self->params->[0] =~ /^s/i);
        }

        $results = {
            CAROL => $result,
        };
    }
    
    return $results;
}

sub compute_carol {

    my ($pph_score, $sift_score) = @_;
    
    my $carol_score;

    if (defined $pph_score) {
        $pph_score = 0.999 if $pph_score == 1;
        $pph_score = 0.0001 if $pph_score == 0;
    }

    if (defined $sift_score) {
        $sift_score = 1 - $sift_score;
        $sift_score = 0.999 if $sift_score == 1;
        $sift_score = 0.0001 if $sift_score == 0;
    }

    if (defined $pph_score && defined $sift_score) {
        
        my $pph_weight  = log(1/(1-$pph_score));
        my $sift_weight = log(1/(1-$sift_score));
       
        # we take -qnorm, because the R script uses qnorm(..., lower.tail = FALSE)

        my $pph_z   = -qnorm($pph_score);
        my $sift_z  = -qnorm($sift_score);
        
        my $numerator   = ($pph_weight * $pph_z) + ($sift_weight * $sift_z);
        my $denominator = sqrt( ($pph_weight ** 2) + ($sift_weight ** 2) );

        # likewise we take 1 - pnorm

        $carol_score = 1 - pnorm($numerator / $denominator);
    }
    elsif (defined $pph_score) {
        $carol_score = $pph_score;
    }
    else {
        $carol_score = $sift_score;
    }

    if (defined $carol_score) {
        my $carol_pred = $carol_score < $CAROL_CUTOFF ? 'Neutral' : 'Deleterious';
        return ($carol_pred, $carol_score);
    }
    else {
        return undef;
    }
}

1;

