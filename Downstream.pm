=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

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

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>
    
=cut

=head1 NAME

 Downstream

=head1 SYNOPSIS

 mv Downstream.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin Downstream

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 predicts the downstream effects of a frameshift variant on the protein
 sequence of a transcript. It provides the predicted downstream protein
 sequence (including any amino acids overlapped by the variant itself),
 and the change in length relative to the reference protein.
 
 Note that changes in splicing are not predicted - only the existing
 translateable (i.e. spliced) sequence is used as a source of
 translation. Any variants with a splice site consequence type are
 ignored.

 If VEP is run in offline mode using the flag --offline, a FASTA file is required.
 See: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#fasta
 Sequence may be incomplete without a FASTA file or database connection.

=cut

package Downstream;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);

    if($self->{config}{offline} && !$self->{config}{fasta}) {
        die("ERROR: cannot function in offline mode without a FASTA file\n");
    }

    return $self;
}

sub version {
    return '2.4';
}

sub feature_types {
    return ['Transcript'];
}

sub variant_feature_types {
    return ['VariationFeature'];
}

sub get_header_info {
    return {
        DownstreamProtein   => "Predicted downstream translation for frameshift mutations",
        ProteinLengthChange => "Predicted change in protein product length",
    };
}

sub run {
    my ($self, $tva) = @_;

    my @SO_terms = map { $_->SO_term } @{$tva->get_all_OverlapConsequences};

    return {} unless grep { $_ eq 'frameshift_variant' } @SO_terms;

    return {} if grep { /splice/ } @SO_terms;

    my $tv = $tva->transcript_variation;
    my $tr = $tv->transcript;

    my $cds_seq = defined($tr->{_variation_effect_feature_cache})
                ? $tr->{_variation_effect_feature_cache}->{translateable_seq}
                : $tr->translateable_seq;

    my ($start, $end) = ($tv->cds_start, $tv->cds_end);

    substr($cds_seq, $start - 1, $end - $start + 1) = $tva->seq_length > 0 ? $tva->feature_seq : '';

    my $low_pos = $start > $end ? $end : $start;
    my $last_complete_codon = $low_pos - ( ( ( $low_pos - 1 ) % 3 ) + 1 );

    my $downstream_seq = substr($cds_seq, $last_complete_codon > 0 ? $last_complete_codon : 0);
    my $three_prime_utr = $tr->three_prime_utr ? $tr->three_prime_utr->seq() : '';

    my $codon_seq = Bio::Seq->new(
      -seq      => $downstream_seq . $three_prime_utr,
      -moltype  => 'dna',
      -alphabet => 'dna'
    );

    my $codon_table;
    if(defined($tr->{_variation_effect_feature_cache})) {
        $codon_table = $tr->{_variation_effect_feature_cache}->{codon_table} || 1;
    }
    else {
        my ($attrib) = @{$tr->slice->get_all_Attributes('codon_table')};
        $codon_table = $attrib ? $attrib->value || 1 : 1;
    }

    my $new_pep = $codon_seq->translate(undef, undef, undef, $codon_table)->seq();
    $new_pep =~ s/\*.*//;

    my $translation = defined($tr->{_variation_effect_feature_cache}->{peptide})
                    ? $tr->{_variation_effect_feature_cache}->{peptide}
                    : $tr->translation->seq;

    my ($pep_start, $pep_end) = ($tv->translation_start, $tv->translation_end);

    my $new_length = ($pep_start < $pep_end ? $pep_start : $pep_end) + length($new_pep);

    return {
        DownstreamProtein   => $new_pep,
        ProteinLengthChange => $new_length - length($translation),
    };
}

1;

