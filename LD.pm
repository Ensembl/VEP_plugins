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

 LD

=head1 SYNOPSIS

 mv LD.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin LD,1000GENOMES:pilot_1_CEU_low_coverage_panel,0.8

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 finds variants in linkage disequilibrium with any overlapping existing 
 variants from the Ensembl variation databases. You can configure the 
 population used to calculate the r2 value, and the r2 cutoff used by 
 passing arguments to the plugin via the VEP command line (separated 
 by commas). This plugin adds a single new entry to the Extra column 
 with a comma-separated list of linked variant IDs and the associated 
 r2 values, e.g.:

  LinkedVariants=rs123:0.879,rs234:0.943

 If no arguments are supplied, the default population used is the CEU
 sample (60 unrelated individuals) from the 1000 Genomes Project pilot 1 low 
 coverage study, and the default r2 cutoff used is 0.8.

 WARNING: Calculating LD is a relatively slow procedure, so this will 
 slow the VEP down considerably when running on large numbers of
 variants. Consider using a filter plugin to limit the analysis to
 'interesting' variants first. You can do this by supplying a filter
 plugin as an argument to the VEP before this one, e.g.:

  perl variant_effect_predictor.pl -i variations.vcf --plugin MyFilter --plugin LD

=cut

package LD;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub version {
    return '2.3';
}

sub feature_types {
    return ['Transcript', 'RegulatoryFeature', 'MotifFeature'];
}

sub get_header_info {

    my $self = shift;
    
    return {
        LinkedVariants => "Variants in LD (r2 >= ".$self->{r2_cutoff}.
            ") with overlapping existing variants from the ".
            $self->{pop}->name." population",
    };
}

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);

    if ($self->config->{cache}) {
        warn "Warning: a connection to the database is required to calculate LD\n";
    }

    my $reg = 'Bio::EnsEMBL::Registry';

    # turn on the check for existing variants

    $self->config->{check_existing} = 1;

    # fetch our population

    my ($pop_name, $r2_cutoff) = @{ $self->params };

    # set some defaults

    $pop_name ||= '1000GENOMES:pilot_1_CEU_low_coverage_panel';

    $r2_cutoff = 0.8 unless defined $r2_cutoff;

    my $pop_adap = $reg->get_adaptor('human', 'variation', 'population')
        || die "Failed to get population adaptor\n";

    my $pop = $pop_adap->fetch_by_name($pop_name)
        || die "Failed to fetch population for name '$pop_name'\n";

    $self->{pop} = $pop;
    $self->{r2_cutoff} = $r2_cutoff;
    
    # prefetch the necessary adaptors
    
    my $ld_adap = $reg->get_adaptor('human', 'variation', 'ldfeaturecontainer')
        || die "Failed to get LD adaptor\n";
    
    my $var_adap = $reg->get_adaptor('human', 'variation', 'variation')
        || die "Failed to get variation adaptor\n";
        
    my $var_feat_adap = $reg->get_adaptor('human', 'variation', 'variationfeature')
        || die "Failed to get variation feature adaptor\n";
     
    $self->{ld_adap} = $ld_adap;
    $self->{var_adap} = $var_adap;
    $self->{var_feat_adap} = $var_feat_adap;

    return $self;
}

sub run {
    my ($self, $vfoa, $line_hash) = @_;

    # fetch the existing variants from the line hash

    my @vars = split ',', $line_hash->{Existing_variation};

    my @linked;

    for my $var (@vars) {
        
        # fetch a variation for each overlapping variant ID

        if (my $v = $self->{var_adap}->fetch_by_name($var)) {

            # and fetch the associated variation features

            for my $vf (@{ $self->{var_feat_adap}->fetch_all_by_Variation($v) }) {

                # we're only interested in variation features that overlap our variant

                if ($vf->slice->name eq $vfoa->variation_feature->slice->name) {

                    # fetch an LD feature container for this variation feature and our preconfigured population

                    if (my $ldfc = $self->{ld_adap}->fetch_by_VariationFeature($vf, $self->{pop})) {
                    
                        # loop over all the linked variants

                        for my $result (@{ $ldfc->get_all_r_square_values }) {
                        
                            # apply our r2 cutoff

                            if ($result->{r2} >= $self->{r2_cutoff}) {

                                my $v1 = $result->{variation1};
                                my $v2 = $result->{variation2};

                                # I'm not sure which of these are the query variant, so just check the names
                                    
                                my $linked = $v1->variation_name eq $var ? $v2 : $v1;
                                
                                push @linked, sprintf("%s:%.3f", $linked->variation_name, $result->{r2});
                            }
                        }
                    }
                }
            }
        }
    }

    # concatenate all our linked variants together

    my $results = join ',', @linked;

    return $results ? {LinkedVariants => $results} : {};
}

1;

