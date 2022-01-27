=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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
 http://ftp.ebi.ac.uk/pub/databases/intact/current/various
 
 2) The genomic location mapped file needs to be tabix indexed. You can 
 do this by following commands -

  a) filter, sort and then zip
  grep -v -e '^$' -e '^[#\-]' mutation_gc_map.txt | sed '1s/.*/#&/' | sort -k1,1 -k2,2n -k3,3n | bgzip > mutation_gc_map.txt.gz
  
  b) create tabix indexed file -
  tabix -s 1 -b 2 -e 3 -f mutation_gc_map.txt.gz

 3) As you have already noticed, tabix utility must be installed in your path to use this plugin.
 
 Options are passed to the plugin as key=value pairs:

 mapping_file			: (mandatory) Path to tabix-indexed genomic location mapped file
 mutation_file			: (mandatory) Path to IntAct data file
 default			: (redundant) Set value to 1 to include feature_type and interaction_ac in the output.
 minimal		  	: Set value to 1 to overwrite default option to include only interaction_ac 
				  in the output
 
 you can also customize output using the following options -
 feature_ac			: Set value to 1 to include Feature AC in the output
 feature_short_label		: Set value to 1 to include Feature short label in the output
 feature_range			: Set value to 1 to include Feature range(s) in the output
 original_sequence		: Set value to 1 to include Original sequence in the output
 resulting_sequence		: Set value to 1 to include Resulting sequence in the output
 feature_type			: Set value to 1 to include Feature type in the output
 feature_annotation		: Set value to 1 to include Feature annotation in the output
 ap_ac				: Set value to 1 to include Affected protein AC in the output
 ap_symbol			: Set value to 1 to include Affected protein symbol in the output
 ap_full_name			: Set value to 1 to include Affected protein full name in the output
 interaction_participants	: Set value to 1 to include Interaction participants in the output
 pmid				: Set value to 1 to include PubMedID in the output
 figure_legend			: Set value to 1 to include Figure legend in the output
 interaction_ac			: (redundant) Set value to 1 to include Interaction AC in the output. Adding
 				  this option redundant as it is already included by default or mimimal options.
 
 See what this options mean - https://www.ebi.ac.uk/intact/download/datasets#mutations
 
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

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();
  
  my $param_hash = $self->params_to_hash();

  $self->add_file($param_hash->{mapping_file});
  
  $self->{interaction_ac} = 1;
  
  $self->{feature_type} = 1 unless defined $param_hash->{minimal};

  $self->{pmid} = 1 if defined $param_hash->{pmid};
  $self->{interaction_participants} = 1 if defined $param_hash->{interaction_participants};   
  $self->{feature_short_label} = 1 if defined $param_hash->{feature_short_label};
  $self->{feature_range} = 1 if defined $param_hash->{feature_range};
  $self->{original_sequence} = 1 if defined $param_hash->{original_sequence};
  $self->{resulting_sequence} = 1 if defined $param_hash->{resulting_sequence};
  $self->{feature_type} = 1 if defined $param_hash->{feature_type};
  $self->{feature_annotation} = 1 if defined $param_hash->{feature_annotation};
  $self->{ap_ac} = 1 if defined $param_hash->{ap_ac};
  $self->{ap_symbol} = 1 if defined $param_hash->{ap_symbol};
  $self->{ap_full_name} = 1 if defined $param_hash->{ap_full_name};
  $self->{figure_legend} = 1 if defined $param_hash->{figure_legend};

  if($self->{config}->{output_format} eq "vcf") {
    $output_vcf = 1;
  }

  if($self->{config}->{output_format} eq "json"){
    $output_json = 1;
  }

  if($self->{config}->{rest}){
    $output_rest = 1;
  }

  return $self;
}

sub feature_types {
  return ['Transcript'];

}

sub get_header_info {
  my $self = shift;
  
  my (%header, %field_des);

  $field_des{"feature_ac"} = "Feature AC - Accession number for that particular mutation feature. ";
  $field_des{"feature_short_label"} = "Feature short label - Human-readable short label summarizing the amino acid changes and their positions (HGVS compliant). ";
  $field_des{"feature_range"} = "Feature range(s) - Position(s) in the protein sequence affected by the mutation. ";
  $field_des{"original_sequence"} = "Original sequence - Wild type amino acid residue(s) affected in one letter code. ";
  $field_des{"resulting_sequence"} = "Resulting sequence - Replacement sequence (or deletion) in one letter code. ";
  $field_des{"feature_type"} = "Feature type - Mutation type following the PSI-MI controlled vocabularies. ";
  $field_des{"feature_annotation"} = "Feature annotation - Specific comments about the feature that can be of interest. ";
  $field_des{"ap_ac"} = "Affected protein AC - Affected protein identifier (preferably UniProtKB accession, if available). ";
  $field_des{"ap_symbol"} = "Affected protein symbol - As given by UniProtKB. ";
  $field_des{"ap_full_name"} = "Affected protein full name - As given by UniProtKB. ";
  $field_des{"interaction_participants"} = "Interaction participants- Identifiers for all participants in the affected interaction along with their species and molecule type between brackets. ";
  $field_des{"pmid"} = "PubMedID - Reference to the publication where the interaction evidence was reported. ";
  $field_des{"figure_legend"} = "Figure legend - Reference to the specific figures in the paper where the interaction evidence was reported. ";
  $field_des{"interaction_ac"} = "Interaction AC - Interaction accession within IntAct databases. ";

  $header{"IntAct"} = "Molecular interaction data from IntAct database. Output may contain multiple interaction data separated by |. Fields in each interaction data are separated by ,. Output field includes :- " unless $output_vcf;

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

  my $param_hash = $self->params_to_hash();
  my $mutation_file = $param_hash->{mutation_file};

  my @matches;

  open my $f, $mutation_file or die("ERROR: couldn't open file $mutation_file\n");
  while(<$f>) {
    my @fields = split /\t/;
    my $ref = $fields[8];
    my (undef, $des) = split /:/, $fields[1];
    
    next unless ($ref and $des);
       
    if( ($des eq $id_des) && ($ref eq $id_ref) ) {
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
  
  # Replace , and | to avoid confusion for interaction paritcipants
  $parsed_data{interaction_participants} =~ s/\|/ and /g;
  $parsed_data{interaction_participants} =~ s/,//g;
   
  return \%parsed_data;
}

# remove duplicate interactions
sub _remove_duplicates {
  my ($self, $intact_matches) = @_;
  
  my @interaction_acs;
  my @feature_types;
  my @interaction_participants;

  my @matches;
  foreach (@$intact_matches) {
    my $parsed_data = $self->_parse_intact_data($_);

    next if $parsed_data->{ap_organism} ne "9606 - Homo sapiens";
    
    my $interaction_ac = $parsed_data->{interaction_ac};
    my $feature_type = $parsed_data->{feature_type};
    my $interaction_participant = $parsed_data->{interaction_participants};
    
    if( grep { $interaction_ac } @interaction_acs &&
        grep { $feature_type } @feature_types &&
        grep { $interaction_participant } @interaction_participants) {

      push @matches, $parsed_data;
    
      push @interaction_acs, $interaction_ac;
      push @feature_types, $feature_type;
      push @interaction_participants, $interaction_participant;
    }
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
          $str = $str ? $str.",".$field_val : $field_val;
        }
      }
    }

    if($output_json || $output_rest){
      my %store_hash = %hash;
      push @arr, \%store_hash;
    }

    unless($output_json || $output_rest || $output_vcf){
      $hash{"IntAct"} = $hash{"IntAct"} ? $hash{"IntAct"}."|".$str : $str; 
    }
  }

  my $filtered_result = $output_json || $output_rest ? {"IntAct" => \@arr} : \%hash;
  
  return $filtered_result;
}

sub run {
  my ($self, $tva) = @_;

  return {} if $self->config->{species} ne "homo_sapiens";

  my $vf = $tva->variation_feature;
  my $tv = $tva->transcript_variation;

  # get reference codon and HGVSp  
  my $ref_codon = $tv->get_reference_TranscriptVariationAllele->codon if $tv;
  my $hgvs = $tva->hgvs_protein;  
  
  return {} unless ($ref_codon and $hgvs);

  my @data =  @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})};
  
  foreach (@data) {
    my ($id_ref, $id_des) = split /:/, $_->{id};
    my (undef, $hgvs_des) = split /:/, $hgvs;
    
    # match variation
    next unless ($id_des eq $hgvs_des and $_->{ref} eq $ref_codon);

    # get matched lines from IntAct data file on hgvs id     
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
