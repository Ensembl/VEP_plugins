=head1 LICENSE
                                                                                                                     
 Copyright (c) 1999-2011 The European Bioinformatics Institute and                                                   
 Genome Research Limited.  All rights reserved.                                                                      
                                                                                                                     
 This software is distributed under a modified Apache license.                                                       
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html                                                               
                                                                                                                     
=head1 CONTACT                                                                                                       

 Will McLaren <wm2@ebi.ac.uk>
    
=cut

=head1 NAME

 Donwstream

=head1 SYNOPSIS

 mv Downstream.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin Downstream

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

=cut

package Downstream;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub version {
    return '2.3';
}

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {
        DownstreamProtein   => "Predicted downstream translation for frameshift mutations",
        ProteinLengthChange => "Predicted change in protein product length",
    };
}

sub run {
    my ($self, $tva) = @_;
    
    my @ocs = @{$tva->get_all_OverlapConsequences};
    
    if(grep {$_->SO_term eq 'frameshift_variant'} @ocs) {
        
        # can't do it for splice sites
        return {} if grep {$_->SO_term =~ /splice/} @ocs;
        
        my $tv = $tva->transcript_variation;
        my $tr = $tv->transcript;
        my $cds_seq = defined($tr->{_variation_effect_feature_cache}) ? $tr->{_variation_effect_feature_cache}->{translateable_seq} : $tr->translateable_seq;
        
        # get the sequence to translate
        my ($low_pos, $high_pos) = sort {$a <=> $b} ($tv->cds_start, $tv->cds_end);
        my $is_insertion         = $tv->cds_start > $tv->cds_end ? 1 : 0;
        my $last_complete_codon  = int($low_pos / 3) * 3;
        my $before_var_seq       = substr $cds_seq, $last_complete_codon, $low_pos - $last_complete_codon - ($is_insertion ? 0 : 1);
        my $after_var_seq        = substr $cds_seq, $high_pos - ($is_insertion ? 1 : 0);
        my $to_translate         = $before_var_seq.$tva->feature_seq.$after_var_seq;
        $to_translate            =~ s/\-//g;
        
        # create a bioperl object
        my $codon_seq = Bio::Seq->new(
          -seq      => $to_translate,
          -moltype  => 'dna',
          -alphabet => 'dna'
        );
        
        # get codon table
        my $codon_table;
        if(defined($tr->{_variation_effect_feature_cache})) {
            $codon_table = $tr->{_variation_effect_feature_cache}->{codon_table} || 1;
        }
        else {
            my ($attrib) = @{$tr->slice->get_all_Attributes('codon_table')};
            $codon_table = $attrib ? $attrib->value || 1 : 1;
        }
        
        # translate
        my $new_pep = $codon_seq->translate(undef, undef, undef, $codon_table)->seq();
        $new_pep =~ s/\*.*//;
        
        # compare lengths
        my $translation = defined($tr->{_variation_effect_feature_cache}) && defined($tr->{_variation_effect_feature_cache}->{peptide}) ? $tr->{_variation_effect_feature_cache}->{peptide} : $tr->translation->seq;
        my $new_length = ($tv->translation_start < $tv->translation_end ? $tv->translation_start : $tv->translation_end) + length($new_pep);
        
        return {
            Downstream          => $new_pep,
            ProteinLengthChange => $new_length - length($translation),
        };
    }

    return {};
}

1;

