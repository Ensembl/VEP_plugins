=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

 Will McLaren <wm2@ebi.ac.uk>
    
=cut

=head1 NAME

 Draw

=head1 SYNOPSIS

 mv SameCodon.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin SameCodon

=head1 DESCRIPTION

 A VEP plugin that reports existing variants that fall in the same codon.

=cut

package SameCodon;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use Bio::EnsEMBL::Variation::Utils::VEP qw(load_dumped_variation_cache);

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub version {
    return '2.5';
}

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {
        SameCodon => "Existing variant IDs that fall in the same codon",
    };
}

sub run {
    my ($self, $tva) = @_;
    
    my $tv = $tva->transcript_variation;
    my $vf = $tv->variation_feature;
    my ($pep_start, $pep_end) = ($tv->translation_start, $tv->translation_end);
    my ($vf_start, $vf_end) = ($vf->start, $vf->end);
    
    return {} unless defined($pep_start) && defined($pep_end);
    
    my $config = $self->{config};
    
    # we need to map the TV start and end coords to the genome
    # needs to be done through the mapper in case the codon spans exons
    my $mapper = $tv->_mapper();
    
    return {} unless defined($mapper);
    
    my @coords = $mapper->pep2genomic($pep_start, $pep_end);
    
    return {} unless scalar @coords;
    return {} if grep {!$_->isa('Bio::EnsEMBL::Mapper::Coordinate')} @coords;
    
    my @results;
    
    # we might get multiple "slices" if the codon that the variant falls in spans exons
    foreach my $coord(@coords) {
        
        my ($slice_start, $slice_end) = ($coord->start, $coord->end);
        
        # using cache?
        if(defined($config->{cache})) {
            
            # spoof region based on cache region size
            my $size = $config->{cache_region_size};
            
            my $s = (int ($vf_start / $size) * $size) + 1;
            my $e = (int ($vf_end / $size) + 1) * $size;
            my $c = $vf->{chr};
            
            my $region = "$s\-$e";
            
            my $vf_cache = $self->{_cache}->{$c}->{$region} || load_dumped_variation_cache($config, $c, $region);
            
            foreach my $pos($slice_start..$slice_end) {
                if(my $existing_vars = $vf_cache->{$c}->{$pos}) {
                    push @results,
                        map {$_->{variation_name}}
                        grep {
                            $_->{variation_name} ne $vf->variation_name &&
                            $_->{failed} <= $config->{failed} &&
                            $_->{start} != $vf_start &&
                            $_->{end} != $vf_end &&
                            scalar $mapper->genomic2cds($_->{start}, $_->{end}, 1) >= 1
                        }
                        @$existing_vars;
                }
            }
            
            $self->{_cache} = {};
            $self->{_cache}->{$c}->{$region} = $vf_cache;
        }
        
        # using DB
        elsif(defined($config->{sa}) && defined($config->{vfa}) && defined($config->{vfa})) {
            
            my $sub_slice = $vf->slice->sub_Slice($slice_start, $slice_end);
            
            push @results,
                map {$_->variation_name}
                grep {
                    $_->variation_name ne $vf->variation_name &&
                    $_->seq_region_start != $vf_start &&
                    $_->seq_region_end != $vf_end &&
                    scalar $mapper->genomic2cds($_->seq_region_start, $_->seq_region_end, 1) >= 1
                }
                @{$sub_slice->get_all_VariationFeatures()};
        }
        
        else {
            return {}
        }
    }
    
    return {} unless scalar @results;
    
    return {
        SameCodon => join ",", @results
    }
}

1;

