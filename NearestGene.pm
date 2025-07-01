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

 NearestGene

=head1 SYNOPSIS

 mv NearestGene.pm ~/.vep/Plugins
 ./vep -i variations.vcf --cache --plugin NearestGene

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 finds the nearest gene(s) to a non-genic variant. More than one gene
 may be reported if the genes overlap the variant or if genes are
 equidistant or if option 'both_directions' is used.

 Various key=value parameters can be altered by passing them to the plugin command:

   limit           : limit the number of genes returned (default: 1)
   range           : initial search range in bp (default: 1000)
   max_range       : maximum search range in bp (default: 50000)
   both_directions : return the nearest genes upstream and downstream of the variant
                     this option overwrites the limit to 1
                     note that the max_range affects the search range in both directions

 Parameters are passed e.g.:

 --plugin NearestGene,limit=3,max_range=50000
 --plugin NearestGene,max_range=50000,both_directions=1

 This plugin requires a database connection. It cannot be run in offline mode i.e. using 
 the --offline flag.

=cut

package NearestGene;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my $char_sep = "|";

my %CONFIG = (
  limit => 1,
  range => 1000,
  max_range => 50000,
);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  
  my $params = $self->params;

  # get output format
  $char_sep = ":" if ($self->{config}->{output_format} eq 'vcf');

  foreach my $param(@$params) {
    my ($key, $val) = split('=', $param);

    die("ERROR: Failed to parse parameter $param\n") unless defined($key) && defined($val);
    if($key eq "both_directions") {
      $self->{both_directions} = 1;
    }
    else {
      $CONFIG{$key} = $val;
    }
  }

  return $self;
}

sub feature_types {
  return ['Intergenic','MotifFeature','RegulatoryFeature'];
}

sub variant_feature_types {
  return ['BaseVariationFeature'];
}

sub get_header_info {
  return {
    NearestGene => "Ensembl identifier of nearest gene"
  };
}

sub run {
  my ($self, $vfoa) = @_;
  
  if ($self->{config}->{offline}) {
      die("ERROR: the plugin NearestGene does not work in --offline mode\n");
  }

  my $vf = $vfoa->base_variation_feature;
  my $loc_string = sprintf("%s:%i-%i", $vf->{chr} || $vf->seq_region_name, $vf->{start}, $vf->{end});
  
  if(!exists($self->{_cache}) || !exists($self->{_cache}->{$loc_string})) {
    $self->{config}->{ga} = $self->{config}->{reg}->get_adaptor($self->{config}->{species}, $self->{config}->{core_type}, 'gene');
    $self->{ga} ||= $self->{config}->{ga};
    die("ERROR: Could not get gene adaptor\n") unless $self->{ga};

    my %opts = map {'-'.$_ => $CONFIG{$_}} keys %CONFIG;
    $opts{-feature} = $vf;

    my @result;

    if($self->{both_directions}) {
      # Overwrite the limit - we want to return only one gene on each direction
      $opts{-limit} = 1;

      # Get upstream genes
      $opts{-upstream} = "upstream";

      my $list_of_genes = $self->{ga}->fetch_all_by_outward_search(%opts);

      for my $gene_result (@{$list_of_genes}){
        my $gene_id = @{$gene_result}[0]->stable_id;
        my $distance = @{$gene_result}[1];
        push(@result, $gene_id.$char_sep.$distance)
      }

      # Get downstream genes
      delete $opts{-upstream};
      $opts{-downstream} = "downstream";

      my $list_of_genes_2 = $self->{ga}->fetch_all_by_outward_search(%opts);

      for my $gene_result (@{$list_of_genes_2}){
        my $gene_id = @{$gene_result}[0]->stable_id;
        my $distance = @{$gene_result}[1];
        push(@result, $gene_id.$char_sep.$distance)
      }
    }
    else {
      # Default behaviour
      @result = map {$_->[0]->stable_id} @{
        $self->{ga}->fetch_all_by_outward_search(%opts)
      };
    }

    $self->{_cache}->{$loc_string} = scalar @result ? join(",", @result) : undef;
  }

  return $self->{_cache}->{$loc_string} ? { NearestGene => $self->{_cache}->{$loc_string} } : {};
}

1;
