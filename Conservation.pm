=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

 Conservation

=head1 SYNOPSIS

 mv Conservation.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin Conservation,GERP_CONSERVATION_SCORE,mammals

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 retrieves a conservation score from the Ensembl Compara databases
 for variant positions. You can specify the method link type and
 species sets as command line options, the default is to fetch GERP
 scores from the EPO 35 way mammalian alignment (please refer to the
 Compara documentation for more details of available analyses). 

 If a variant affects multiple nucleotides the average score for the
 position will be returned, and for insertions the average score of
 the 2 flanking bases will be returned.

 The plugin requires the ensembl-compara API module to be installed;
 see http://www.ensembl.org/info/docs/api/index.html

=cut

package Conservation;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Slice;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub version {
    return '2.3';
}

sub feature_types {
    return ['Feature','Intergenic'];
}

sub get_header_info {

    my $self = shift;
    
    return {
        Conservation => "The conservation score for this site (method_link_type=\"".
            $self->{method_link_type}."\", species_set=\"".$self->{species_set}."\")",
    };
}

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);

    my $params = $self->params;

    # REST API passes 1 as first param
    shift @$params if $params->[0] && $params->[0] eq '1';

    $self->{method_link_type} = $params->[0] || 'GERP_CONSERVATION_SCORE';
    $self->{species_set}  = $params->[1] || 'mammals';

    my $config = $self->{config};
    my $reg = 'Bio::EnsEMBL::Registry';

    # reconnect to DB without species param
    if($config->{host}) {
        $reg->load_registry_from_db(
            -host       => $config->{host},
            -user       => $config->{user},
            -pass       => $config->{password},
            -port       => $config->{port},
            -db_version => $config->{db_version},
            -no_cache   => $config->{no_slice_cache},
        );
    }

    my $mlss_adap = $config->{mlssa} ||= $reg->get_adaptor('Multi', 'compara', 'MethodLinkSpeciesSet')
        or die "Failed to connect to compara database\n";

    $self->{mlss} = $mlss_adap->fetch_by_method_link_type_species_set_name($self->{method_link_type}, $self->{species_set})
        or die "Failed to fetch MLSS for ".$self->{method_link_type}." and ".$self->{species_set}."\n";

    $self->{cs_adap} = $config->{cosa} ||= $reg->get_adaptor('Multi', 'compara', 'ConservationScore')
        or die "Failed to fetch conservation adaptor\n";

    return $self;
}

sub run {
    my ($self, $bvfoa) = @_;

    my $bvf = $bvfoa->base_variation_feature;

    # we cache the score on the BaseVariationFeature so we don't have to
    # fetch it multiple times if this variant overlaps multiple Features

    unless (exists $bvf->{_conservation_score}) {

        my $slice;

        my $true_snp = 0;

        if ($bvf->{end} >= $bvf->{start}) {

            if ($bvf->{start} == $bvf->{end}) {

                # work around a bug in the compara API that means you can't fetch 
                # conservation scores for 1bp slices by creating a 2bp slice for
                # SNPs and then ignoring the score returned for the second position

                my $s = $bvf->slice;

                $slice = Bio::EnsEMBL::Slice->new(
                    -seq_region_name   => $s->seq_region_name,
                    -seq_region_length => $s->seq_region_length,
                    -coord_system      => $s->coord_system,
                    -start             => $bvf->{start},
                    -end               => $bvf->{end} + 1,
                    -strand            => $bvf->{strand},
                    -adaptor           => $s->adaptor
                );
                
                $true_snp = 1;
            }
            else {

                # otherwise, just get a slice that covers our variant feature

                $slice = $bvf->feature_Slice;
            }
        }
        else {

            # this is an insertion, we return the average score of the flanking 
            # bases, so we create a 2bp slice around the insertion site

            my $s = $bvf->slice;

            $slice = Bio::EnsEMBL::Slice->new(
                -seq_region_name   => $s->seq_region_name,
                -seq_region_length => $s->seq_region_length,
                -coord_system      => $s->coord_system,
                -start             => $bvf->{end},
                -end               => $bvf->{start},
                -strand            => $bvf->{strand},
                -adaptor           => $s->adaptor
            );
        }

        my $scores = $self->{cs_adap}->fetch_all_by_MethodLinkSpeciesSet_Slice(
            $self->{mlss},                      # our MLSS for the conservation metric and the set of species
            $slice,                             # our slice
            ($slice->end - $slice->start + 1),  # the number of scores we want back (one for each base)
        );

        if (@$scores > 0) {

            # we use the simple average of the diff_scores as the overall score
            
            pop @$scores if $true_snp; # get rid of our spurious second score for SNPs
            
            if (@$scores > 0) {

                my $tot_score = 0;
    
                $tot_score += $_->diff_score for @$scores;
    
                $tot_score /= @$scores;
                
                $bvf->{_conservation_score} = sprintf "%.3f", $tot_score;
            }
            else {
                $bvf->{_conservation_score} = undef;
            }
        }
        else {
            $bvf->{_conservation_score} = undef;
        }
    }

    if (defined $bvf->{_conservation_score}) {
        return {
            Conservation => $bvf->{_conservation_score}
        };
    }
    else {
        return {};
    }
}

1;

