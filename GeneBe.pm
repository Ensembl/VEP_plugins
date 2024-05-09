=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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
 Piotr Stawinski <genebeambio@gmail.com>

=cut

=head1 NAME

 GeneBe - ACMG automatic criteria assignment

=head1 SYNOPSIS

 mv GeneBe.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin GeneBe

=head1 DESCRIPTION

 A user-contributed VEP plugin that retrieves automatic ACMG variant classification data from
 https://genebe.net/

 Please cite the GeneBe publication alongside the VEP if you use this resource:
 https://onlinelibrary.wiley.com/doi/10.1111/cge.14516 .

 Please be advised that the GeneBe API is freely accessible for academic purposes only, with a limited
 number of queries per day, albeit at a high threshold. Kindly utilize this resource judiciously
 to ensure its availability for others. For further information, please visit https://genebe.net/about/api.

=cut

package GeneBe;

use strict;
use warnings;
use JSON;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub version {
    return '1.0';
}

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {
        GeneBe_ACMG_score => "ACMG score automatically assigned by GeneBe algorithm",
        GeneBe_ACMG_classification => "ACMG classification automatically assigned by GeneBe algorithm",
        GeneBe_ACMG_criteria => "ACMG criteria applied by GeneBe algorithm",
    };
}

sub run {
    my ($self, $tva) = @_;

    $self->{has_cache} = 1;

    my %hg_assembly = ( #'grch37' => 'hg19',
        'grch38' => 'hg38' );
    my $assembly = lc($self->{config}->{assembly});

    # only works on human
    die("ERROR: GeneBe ACMG plugin works only on human data") unless $self->{config}->{species} =~ /human|homo/i;

    # only work on 1 human assembly right now
    die("ERROR: GeneBe ACMG plugin works only on human assembly hg38 ".join(' and ', keys(%hg_assembly))) unless $hg_assembly{$assembly};



    # get the VF object
    my $vf = $tva->variation_feature;
    my $transcript = $tva->transcript;
    my $end = $vf->{end};
    my $start = $vf->{start};

    my $ref_allele = $vf->ref_allele_string;
    my $alt_alleles = $tva->base_variation_feature->alt_alleles;
    my $allele_number = $tva->{allele_number};
    my $alt_allele = @$alt_alleles[$allele_number - 1];

    return {} unless defined $vf;

    my $feature = $vf;

    # Convert LRG to chromosome system
    if ( $feature->coord_system_name() eq "lrg" ) {
        $feature = $vf->transform('chromosome');
        unless ( $feature ) {
            warn "Region $vf->{chr}:$vf->{start}_$vf->{end} not defined in " .
                 "chromosome coordinate system\n";
            return {};
        }
    }
    my $chr = $feature->slice->seq_region_name();
    $chr =~ s/^chr//;

    my $locus = sprintf 'chr%s-%s-%s-%s', $chr, $start, $ref_allele, $alt_allele;

    # check the cache
    if(!exists($self->{genebe_cache}->{$locus})) {

        # Define the base URL of your API
        my $base_url = 'https://api.genebe.net/cloud/api-public/v1/variant';

        # Construct the URL with the variables
        my $api_url = sprintf('%s?chr=%s&pos=%d&ref=%s&alt=%s&genome=hg38',
                      $base_url, $chr, $start, $ref_allele, $alt_allele);

        # Execute curl command to make the HTTP request
        my $curl_command = "curl --netrc -s -A GeneBe_VEP_plugin \"$api_url\"";
        my $curl_output = `$curl_command`;

        if($? == 0) {

            # Decode the JSON response
            my $json_response = decode_json($curl_output);

            # Check if 'variants' array exists and has at least one element
            if (exists $json_response->{variants} && @{$json_response->{variants}}) {
                # Extract data from the first element
                my $variant_data = $json_response->{variants}[0];

                # Extract required values
                my $acmg_score = $variant_data->{acmg_score} // '';
                my $acmg_classification = $variant_data->{acmg_classification} // '';
                my $acmg_criteria = $variant_data->{acmg_criteria} // '.';
                $acmg_criteria = [ split /,/, $acmg_criteria ] if defined $acmg_criteria;

                $self->{genebe_cache}->{$locus} = [$acmg_score, $acmg_classification, $acmg_criteria];

                return {GeneBe_ACMG_score => $acmg_score, GeneBe_ACMG_classification => $acmg_classification, GeneBe_ACMG_criteria => $acmg_criteria };
            }
        }
    }
    else {
        my $cached_data = $self->{genebe_cache}->{$locus};
        my ($acmg_score, $acmg_classification, $acmg_criteria) = @$cached_data;
        return {GeneBe_ACMG_score => $acmg_score, GeneBe_ACMG_classification => $acmg_classification, GeneBe_ACMG_criteria => $acmg_criteria };
    }

     return {};
}

1;

