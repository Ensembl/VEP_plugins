=head1 LICENSE

GNU General Public License, version 3 (GPL-3.0)

=head1 CONTACT

 Questions may also be sent to <duarte.molha@ogt.co.uk>.

=cut

=head1 NAME

qual_score_and_genotype

=head1 SYNOPSIS

    mv qual_score_and_genotype.pm ~/.vep/Plugins
    perl variant_effect_predictor.pl -i variations.vcf --plugin qual_score_and_genotype
 
 
=head1 DESCRIPTION

This VEP plugin retrieves the genotype, the format genotype and the score fields and writes each argument in key-value pairs in the Extra column of the VEP.
This is usefull if your outputing a tab delimited annotation file from VEP but you are interested in
retaining these values from the original VCF input.
in case your output annotation is a VCF file then this plugin is a bit redundant and should probably be avoided. 

=cut

package qual_score_and_genotype;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

use strict;
use warnings;

sub version {
    return '2.6';
}


sub get_header_info {
    return {
        "quality_score" => "Quality score from VCF input Field",
        "AD" => "Allelic depths for the ref and alt alleles in the order listed",
        "DP" => "Read Depth (only filtered reads used for calling)",
        "GQ" => "Genotype Quality",
        "GT" => "Genotype",
        "PL" => "Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic",
    };
}

sub feature_types {
    return ['Feature', 'Intergenic'];
}


sub run {
    my $self = shift;
    my $vf = shift;
    my $line_hash = shift;
    
    my $config = $self->{config};
    my $ind_cols = $config->{ind_cols};
    my $line = $vf->{base_variation_feature_overlap}->{base_variation_feature}->{_line};
    my $individual = $vf->{base_variation_feature_overlap}->{base_variation_feature}->{individual};

    my @split_line = split /[\s\t]+/, $line;
    my $qual_score = $split_line[5];
    my @gt_format  = split /:/, $split_line[8];
    my @gt_data    = split /:/, $split_line[$ind_cols->{$individual}];
    my $results = {map { shift @gt_format => $_ } @gt_data};
    $results->{"quality_score"} = $qual_score;
    return $results;
}

1;
