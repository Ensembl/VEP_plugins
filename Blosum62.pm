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

 Blosum62

=head1 SYNOPSIS

 mv Blosum62.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin Blosum62

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 looks up the BLOSUM 62 substitution matrix score for the reference
 and alternative amino acids predicted for a missense mutation. It adds
 one new entry to the VEP's Extra column, BLOSUM62 which is the 
 associated score. 

=cut

package Blosum62;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my @BLOSUM_62 = qw(
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2 
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
);

my @AAs = qw(A R N D C Q E G H I L K M F P S T W Y V);

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);

    # construct a hash representing the matrix for quick lookups

    my $num = @AAs;

    for (my $i = 0; $i < $num; $i++) {
        for (my $j = 0; $j < $num; $j++) {
            $self->{matrix}->{$AAs[$i]}->{$AAs[$j]} = $BLOSUM_62[($i * $num) + $j];
        }
    }

    return $self;
}

sub version {
    return '2.3';
}

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {
        BLOSUM62 => "BLOSUM62 substitution score for the reference and alternative amino acids",
    };
}

sub run {
    my ($self, $tva) = @_;

    if ($tva->pep_allele_string && $tva->pep_allele_string =~ /^([A-Z])\/([A-Z])$/) {
        
        my $score = $self->{matrix}->{$1}->{$2};
        
        if (defined $score) {
            return {
                BLOSUM62 => $score
            };
        }
    }

    return {};
}

1;

