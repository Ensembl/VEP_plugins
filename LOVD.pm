=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and                                                   
 Genome Research Limited.  All rights reserved.                                                                      

 This software is distributed under a modified Apache license.                                                       
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html                                                               

=head1 CONTACT                                                                                                       

 Will McLaren <wm2@ebi.ac.uk>
    
=cut

=head1 NAME

 Draw

=head1 SYNOPSIS

 mv LOVD.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin LOVD

=head1 DESCRIPTION

 A VEP plugin that retrieves LOVD variation data from http://www.lovd.nl/.
 
 Please be aware that LOVD is a public resource of curated variants, therefore
 please respect this resource and avoid intensive querying of their databases
 using this plugin, as it will impact the availability of this resource for
 others.

=cut

package LOVD;

use strict;
use warnings;
use LWP::UserAgent;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub version {
    return '2.5';
}

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {
        LOVD => "LOVD variant ID",
    };
}

sub run {
    my ($self, $tva) = @_;
    
    $self->{has_cache} = 1;
    
    # only works on human
    die("ERROR: LOVD plugin works only on human data") unless $self->{config}->{species} =~ /human|homo/i;
    
    # get the VF object
    my $vf = $tva->variation_feature;
    return {} unless defined $vf;
    
    # set up a LWP UserAgent
    my $ua = LWP::UserAgent->new;
    $ua->env_proxy;
    
    my $locus = sprintf('chr%s:%s_%s', $vf->{chr}, $vf->{start}, $vf->{end});
    
    my $data;
    
    # check the cache
    if(!exists($self->{lovd_cache}->{$locus})) {
        
        # construct a LOVD URL
        my $url = 'http://www.lovd.nl/search.php?build=hg19&position='.$locus;
        
        # get the accession (only need the head to get the redirect URL that contains the accession)
        my $response = $ua->get($url);
        
        if($response->is_success) {
            
            # parse the data into a hash
            for(grep {$_ !~ /hg_build/} split /\cJ/, $response->decoded_content) {
                s/\"//g;
                
                my ($build, $pos, $gene, $acc, $dna, $id, $url) = split /\t/;
                
                $data->{$id} = {
                    gene => $gene,
                    acc  => $acc,
                    dna  => $dna
                };
            }
            
            $self->{lovd_cache}->{$locus} = $data;
        }
    }
    else {
        $data = $self->{lovd_cache}->{$locus};
    }
    
    return {} unless scalar keys %$data;
    
    return {
        LOVD => (join ",", keys %$data)
    };
}

1;

