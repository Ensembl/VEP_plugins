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

