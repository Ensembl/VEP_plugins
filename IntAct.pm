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

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>
    
=cut

=head1 NAME

 IntAct

=head1 SYNOPSIS

 mv IntAct.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin IntAct,mutation_file=/FULL_PATH_TO_IntAct_FILE/mutations.tsv,mapping_file=/FULL_PATH_TO_IntAct_FILE/mutation_gc_map.txt.gz
 ./vep -i variations.vcf --plugin IntAct,mutation_file=/FULL_PATH_TO_IntAct_FILE/mutations.tsv,mapping_file=/FULL_PATH_TO_IntAct_FILE/mutation_gc_map.txt.gz,minimal=1

=head1 DESCRIPTION

 A VEP plugin that retrieves molecular interaction data for variants as reprted by IntAct database.
 
 Please cite the IntAct publication alongside the VEP if you use this resource:
 https://pubmed.ncbi.nlm.nih.gov/24234451/
 
 Pre-requisites:
 
 1) IntAct files can be downloaded from -
 https://ftp.ebi.ac.uk/pub/databases/intact/current/various
 
 2) The genomic location mapped file needs to be tabix indexed. You can 
 do this by following commands -

  a) filter, sort and then zip
  grep -v -e '^$' -e '^[#\-]' mutation_gc_map.txt | sed '1s/.*/#&/' | awk -F "\t" 'BEGIN { OFS="\t"} {if ($2 > $3) {a=$2; $2=$3; $3=a}; print $0 }' | sort -k1,1 -k2,2n -k3,3n | bgzip > mutation_gc_map.txt.gz
  
  b) create tabix indexed file -
  tabix -s 1 -b 2 -e 3 -f mutation_gc_map.txt.gz

 3) As you have already noticed, tabix utility must be installed in your path to use this plugin.
 
 Options are passed to the plugin as key=value pairs:

 mapping_file     : (mandatory) Path to tabix-indexed genomic location mapped file
 mutation_file    : (mandatory) Path to IntAct data file
 
 By default the output will always contain feature_type and interaction_ac from the IntAct data file. You can also add more fields using the following key=value options -
   feature_ac               : Set value to 1 to include Feature AC in the output
   feature_short_label      : Set value to 1 to include Feature short label in the output
   feature_annotation       : Set value to 1 to include Feature annotation in the output
   ap_ac                    : Set value to 1 to include Affected protein AC in the output
   interaction_participants : Set value to 1 to include Interaction participants in the output
   pmid                     : Set value to 1 to include PubMedID in the output

 There are also two other key=value options for customizing the output - 
   all     : Set value to 1 to include all the fields
   minimal : Set value to 1 to overwrite default behavior and include only interaction_ac 
				  in the output by default

 See what these options mean - https://www.ebi.ac.uk/intact/download/datasets#mutations
 
 Note that, interaction accession can be used to link to full details on the interaction website. For example, 
 where the VEP output reports an interaction_ac of EBI-12501485, the URL would be : 

   https://www.ebi.ac.uk/intact/details/interaction/EBI-12501485

=cut

package IntAct;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my $output_vcf;
my $output_json;
my $output_rest;

# hashref containing all fields that exists in IntAct data file 
my $valid_fields = {
  0 => "feature_ac",
  1 => "feature_short_label",
  2 => "feature_range",
  3 => "original_sequence",
  4 => "resulting_sequence",
  5 => "feature_type",
  6 => "feature_annotation",
  7 => "ap_ac",
  8 => "ap_symbol",
  9 => "ap_full_name",
  10 => "ap_organism",
  11 => "interaction_participants",
  12 => "pmid",
  13 => "figure_legend",
  14 => "interaction_ac"
};

# taxomy id for supported species
my $taxonomy_lookup = {
  "homo_sapiens" => 9606,
  "mus_musculus" => 10090,
  "rattus_norvegicus" => 10116,
  "gallus_gallus_gca000002315v5" => 9031,
  "saccharomyces_cerevisiae" => 559292,
  "arabidopsis_thaliana" => 3702,
  "drosophila_melanogaster" => 7227
};

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();
  
  my $param_hash = $self->params_to_hash();

  die "ERROR: mutation_file is not specified which is a mandatory parameter\n" unless defined $param_hash->{mutation_file};
  die "ERROR: mapping_file is not specified which is a mandatory parameter\n" unless defined $param_hash->{mapping_file};

  $self->add_file($param_hash->{mapping_file});
  
  $self->{interaction_ac} = 1;
  
  $self->{feature_type} = 1 unless defined $param_hash->{minimal};

  $self->{feature_ac} = 1 if $param_hash->{feature_ac} || $param_hash->{all};   
  $self->{feature_short_label} = 1 if $param_hash->{feature_short_label} || $param_hash->{all};
  $self->{feature_type} = 1 if $param_hash->{all};
  $self->{feature_annotation} = 1 if $param_hash->{feature_annotation} || $param_hash->{all};
  $self->{ap_ac} = 1 if $param_hash->{ap_ac} || $param_hash->{all};
  $self->{interaction_participants} = 1 if $param_hash->{interaction_participants} || $param_hash->{all};
  $self->{pmid} = 1 if $param_hash->{pmid} || $param_hash->{all};

  if($self->{config}->{output_format} eq "vcf") {
    $output_vcf = 1;
  }

  if($self->{config}->{output_format} eq "json"){
    $output_json = 1;
  }

  if($self->{config}->{rest}){
    $output_rest = 1;
  }

  my $species = $self->config->{species};
  $self->{species} = $species;

  # get species taxonomy id using lookup
  my $species_tax_id = $taxonomy_lookup->{$species};
  $self->{species_tax_id} = $species_tax_id;

  die "ERROR: could not get species taxonomy id\n" unless defined $species_tax_id;

  return $self;
}

sub feature_types {
  return ['Transcript'];

}

sub get_header_info {
  my $self = shift;
  
  my (%header, %field_des);

  $field_des{"feature_ac"} = "Feature AC - Accession number for that particular mutation feature; ";
  $field_des{"feature_short_label"} = "Feature short label - Human-readable short label summarizing the amino acid changes and their positions (HGVS compliant); ";
  $field_des{"feature_type"} = "Feature type - Mutation type following the PSI-MI controlled vocabularies; ";
  $field_des{"feature_annotation"} = "Feature annotation - Specific comments about the feature that can be of interest; ";
  $field_des{"ap_ac"} = "Affected protein AC - Affected protein identifier (preferably UniProtKB accession, if available); ";
  $field_des{"interaction_participants"} = "Interaction participants- Identifiers for all participants in the affected interaction along with their species and molecule type between brackets; ";
  $field_des{"pmid"} = "PubMedID - Reference to the publication where the interaction evidence was reported; ";
  $field_des{"interaction_ac"} = "Interaction AC - Interaction accession within IntAct databases; ";

  $header{"IntAct"} = "Molecular interaction data from IntAct database. Output may contain multiple interaction data separated by ,. Fields in each interaction data are separated by |. Output field includes: " unless $output_vcf;

  my $i = 0;
  my $total_fields = scalar keys %$valid_fields;

  while ($i < $total_fields){
    my $field = $valid_fields->{$i++};
    next unless ($self->{$field} and $field ne "ap_organism");
 
    if($output_vcf){
      $header{"IntAct_".$field} = $field_des{$field};
    }
    else{
      $header{"IntAct"} .= $field_des{$field};
    }
  }
   
  return \%header;
}

# match lines of IntAct data file using HGVS id
sub _match_id {
  my ($self, $id_ref, $id_des) = @_;

  # remove any 'p.' prefix if exists
  $id_des =~ s/p\.//;

  my $param_hash = $self->params_to_hash();
  my $mutation_file = $param_hash->{mutation_file};

  my @matches;

  open my $f, $mutation_file or die("ERROR: couldn't open file $mutation_file\n");
  while(<$f>) {
    my @fields = split /\t/;
    my $ref = $fields[8];
    my (undef, $des) = split /:/, $fields[1];
    
    next unless ($ref and $des);
    
    # remove any 'p.' prefix if exists
    $des =~ s/p\.//;

    if( ($des =~ /$id_des/) && ($ref eq $id_ref) ) {
        push @matches, $_;
    }
  }

  return \@matches;
}

# parse the extracted lines of IntAct data file
sub _parse_intact_data {
  my ($self, $data) = @_;

  my $i = 0;
  my %parsed_data = map { $valid_fields->{$i++} => $_ } split /\t/, $data
    or die("ERROR: cannot parse intact data.\n");
  
  # process complex field for readability
  if ( defined $parsed_data{interaction_participants} && $parsed_data{interaction_participants} =~ /[\(\)\|]/){
    $parsed_data{interaction_participants} = join ' and ', map { s/\(.*//g; $_ } split /\|/, $parsed_data{interaction_participants};
  }

  if ( defined $parsed_data{feature_type} && $parsed_data{feature_type} =~ /[\(\)]/ ){
    $parsed_data{feature_type} =~ s/\([^\)]+\)//g;
  }  
 
  return \%parsed_data;
}

# check if the data has correct species
sub _match_species{
  my ($self, $ap_organism) = @_;

  if( defined $self->{species_tax_id} ){
    # get taxonomy id from the ap_organism field
    my $ap_organism_tax_id = (split /-/, $ap_organism)[0];
    $ap_organism_tax_id =~ s/^\s+|\s+$//;

    return 1 if $self->{species_tax_id} eq $ap_organism_tax_id;
  }
  
  return 0;
}

# remove duplicate interactions
sub _remove_duplicates {
  my ($self, $intact_matches) = @_;

  my @matches;
  
  my $taken = {};
  foreach (@$intact_matches) {
    my $is_unique = 0;

    my $parsed_data = $self->_parse_intact_data($_);

    next unless $self->_match_species($parsed_data->{ap_organism});

    foreach (keys %$parsed_data) {
      next unless defined $self->{$_};

      my ($field, $field_val) = ($_, $parsed_data->{$_});
      chomp $field_val;

      $taken->{$field} = [] unless defined $taken->{$field};

      unless( grep { /\Q$field_val/ } @{ $taken->{$field} } ) {
        push @{ $taken->{$field} }, $field_val;
        $is_unique = 1;
      }
    }

    push @matches, $parsed_data if $is_unique;
  }
  
  return \@matches;
}

# filter out field according to user input
sub _filter_fields {
  my ($self, $uniq_matches) = @_;
  
  my (@arr, %hash);
  
  foreach (@$uniq_matches) {
    my $str;

    my $j = 0;
    while (defined $valid_fields->{$j}) {
      my $field = $valid_fields->{$j++};
      my $field_val = $_->{$field};

      chomp $field_val;

      if(defined $self->{$field}) {
        if($output_json || $output_rest){
          $hash{$field} = $field_val;
        }
        elsif($output_vcf){
          $hash{"IntAct_".$field} = $hash{"IntAct_".$field} ? $hash{"IntAct_".$field}.",".$field_val : $field_val;
        }
        else{
          $str = $str ? $str."|".$field_val : $field_val;
        }
      }
    }

    if($output_json || $output_rest){
      my %store_hash = %hash;
      push @arr, \%store_hash;
    }

    unless($output_json || $output_rest || $output_vcf){
      $hash{"IntAct"} = $hash{"IntAct"} ? $hash{"IntAct"}.",".$str : $str; 
    }
  }

  my $filtered_result = $output_json || $output_rest ? {"IntAct" => \@arr} : \%hash;
  
  return $filtered_result;
}

sub run {
  my ($self, $tva) = @_;

  # skip if the species is not supported
  next unless exists $taxonomy_lookup->{ $self->{species} };

  my $vf = $tva->variation_feature;
  my $tv = $tva->transcript_variation;

  # get HGVSp for the variant in question through variation API  
  my $hgvs = $tva->hgvs_protein;
  return {} unless $hgvs;

  my @data =  @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})};
  
  foreach (@data) {
    my ($id_ref, $id_des) = split /:/, $_->{id};
    my (undef, $hgvs_des) = split /:/, $hgvs;
    
    # match HGVSp description to make sure we get the correct variant from data file
    next unless $id_des eq $hgvs_des;

    # get matched lines from IntAct data file on the hgvs like description 
    my $intact_matches = $self->_match_id($id_ref, $id_des);

    # keep only the unique interaction data from IntAct
    my $uniq_matches = $self->_remove_duplicates($intact_matches) if $intact_matches;

    # filter fields according to user defined parameters
    my $results = $self->_filter_fields($uniq_matches) if $uniq_matches;

    return $results if $results;
  }

  return {};
}

sub parse_data {
  my ($self, $line) = @_;
  my ($c, $s, $e, $ref, undef, undef, $id) = split /\t/, $line;
  
  return {
    chr => $c,
    start => $s,
    end => $e,
    ref => $ref,
    id => $id
  };
}

sub get_start { 
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
