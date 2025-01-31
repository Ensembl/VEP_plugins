=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

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

 TSSDistance

=head1 SYNOPSIS

 mv TSSDistance.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin TSSDistance

# Get both up and downstream distances:
 ./vep -i variations.vcf --plugin TSSDistance,both_direction=1

=head1 DESCRIPTION

 A VEP plugin that calculates the distance from the transcription
 start site for upstream variants. Or variants in both directions
 if parameter `both_direction=1` is provided.

=cut

package TSSDistance;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);

    # Get directionality flag:
    my $param_hash = $self->params_to_hash();

    # If both_direction parameter is set, we will calculate the distance from the TSS in both directions, otherwise only upstream distance will be calculated:
    $self->{both_direction} = defined($param_hash->{both_direction}) ? 1 : 0;
    return $self;
}

sub get_header_info {
    return {
        TSSDistance => "Distance from the transcription start site"
    };
}

sub feature_types {
    return ['Transcript'];
}

sub variant_feature_types {
    return ['BaseVariationFeature'];
}

sub run {
    my ($self, $tva) = @_;

    my $t = $tva->transcript;
    my $vf = $tva->base_variation_feature;

    my $dist;

    if ($t->strand == 1) {
        $dist = $t->start - $vf->end;
    }
    else {
        $dist = $vf->start - $t->end;
    }

    # Return upstream distance:
    if ($dist > 0) {
        return {
            TSSDistance => $dist,
        }
    }

    # Downstream distance only returned if both_direction flag is set to 1:
    if ($self->{both_direction} == 1){
        return {
            TSSDistance => $dist,
        }
    }

    # Otherwise return empty hash:
    return {};
}

1;
