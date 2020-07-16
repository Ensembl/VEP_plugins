=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

 neXtProt

=head1 SYNOPSIS

 mv neXtProt.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin neXtProt
 ./vep -i variations.vcf --plugin neXtProt,max_set=1

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 retrives data for missense variants from neXtProt, which is a comprehensive 
 human-centric discovery platform that offers integration of and navigation 
 through protein-related data (https://www.nextprot.org/).

 Please cite the neXtProt publication alongside the VEP if you use this resource:
 https://doi.org/10.1093/nar/gkz995
 
 Running options:
 (Default) the data retrieved by default is the MatureProtein, NucleotidePhosphateBindingRegion,
 Variant, Domain, MiscellaneousRegion and InteractingRegion.
 The plugin can also be run with other options to retrieve other data than the default.
 
 Options are passed to the plugin as key=value pairs:
 max_set        : Set value to 1 to return a bigger set of protein-related data 
                  (includes the default data)

 config_data    : The set of data can be configured by the user. Use file 'neXtProt_headers.txt'
                  to check which data (labels) are available.

 url            : Set value to 1 to include the URL to link to the neXtProt entry

 all_labels     : Set value to 1 to include all data even if data was not found for all labels

 * note: 'max_set' and 'config_data' cannot be used in simultaneously.


 Output:
  By default, the plugin only returns data that is available. Example (default behaviour):
  neXtProt_Domain=396,583,DH;neXtProt_MatureProtein=1,1344,Rho guanine nucleotide exchange factor 10

  The option 'all_labels' includes all data, same example as above:
  neXtProt_Domain=396,583,DH;neXtProt_MatureProtein=1,1344,Rho guanine nucleotide exchange factor 10;
  neXtProt_InteractingRegion=-;neXtProt_NucleotidePhosphateBindingRegion=-;neXtProt_Variant=-;
  neXtProt_MiscellaneousRegion=-;


 The plugin can then be run as default:
 ./vep -i variations.vcf --plugin neXtProt

 or to return only the data specified by the user:
 ./vep -i variations.vcf --plugin neXtProt,config_data='Domain&InteractingRegion'


=cut

package neXtProt;

use strict;
use warnings;
use JSON::XS;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my $default_output = {
  'neXtProt_MatureProtein' => 'Extent of an active peptide or a polypetide chain in the mature protein',
  'neXtProt_NucleotidePhosphateBindingRegion' => 'Nucleotide phosphate binding region',
  'neXtProt_Variant' => 'Natural variant of the protein',
  'neXtProt_Domain' => 'Position and type of each modular protein domain',
  'neXtProt_MiscellaneousRegion' => 'Region of interest in the sequence',
  'neXtProt_InteractingRegion' => 'Region interacting with another macromolecule'
};

my $max_set_output = {
  'neXtProt_TopologicalDomain' => 'Location of non-membrane regions of membrane-spanning proteins',
  'neXtProt_SequenceConflict' => 'Sequence discrepancies of unknown origin',
  'neXtProt_TransmembraneRegion' => 'Extent of a membrane-spanning region',
  'neXtProt_CompositionallyBiasedRegion' => 'Region of compositional bias in the protein',
  'neXtProt_ModifiedResidue' => 'Modified residues',
  'neXtProt_Repeat' => 'Positions of repeated sequence motifs or repeated domains',
  'neXtProt_Mutagenesis' => 'Site which has been experimentally altered by mutagenesis',
  'neXtProt_DisulfideBond' => 'Cysteine residues participating in disulfide bonds',
  'neXtProt_GlycosylationSite' => 'Covalently attached glycan group(s)',
  'neXtProt_ZincFingerRegion' => 'Position(s) and type(s) of zinc fingers within the protein',
  'neXtProt_DnaBindingRegion' => 'Position and type of a DNA-binding domain',
  'neXtProt_BindingSite' => 'Binding site for any chemical group (co-enzyme, prosthetic group, etc.)',
  'neXtProt_MatureProtein' => 'Extent of an active peptide or a polypetide chain in the mature protein',
  'neXtProt_NucleotidePhosphateBindingRegion' => 'Nucleotide phosphate binding region',
  'neXtProt_Variant' => 'Natural variant of the protein',
  'neXtProt_Domain' => 'Position and type of each modular protein domain',
  'neXtProt_MiscellaneousRegion' => 'Region of interest in the sequence',
  'neXtProt_InteractingRegion' => 'Region interacting with another macromolecule'
};

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $param_hash = $self->params_to_hash();

  if(defined($param_hash->{max_set}) && defined($param_hash->{config_data})) {
    die "ERROR: Can't use max_set and config_data simultaneously!\n";
  }

  # Return the isoform URL to the neXtProt web page
  if(defined($param_hash->{url})) {
    $self->{url} = $param_hash->{url};

    if(defined($param_hash->{max_set})) {
      $max_set_output->{'neXtProt_url'} = 'neXtProt URL';
    }
    elsif(defined($param_hash->{config_data})) {
      $self->{config_data_hash}->{'neXtProt_url'} = 'neXtProt URL';
    }
    else {
      $default_output->{'neXtProt_url'} = 'neXtProt URL';
    }
  }

  if(defined($param_hash->{max_set})) {
    $self->{max_set} = $param_hash->{max_set};
  }

  if(defined($param_hash->{config_data})) {
    $self->{config_data} = $param_hash->{config_data};
    $self->build_data_hash();
  }

  if(defined($param_hash->{all_labels})) {
    $self->{all_labels} = $param_hash->{all_labels};
  }

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  my %header;

  if($self->{max_set}) {
    foreach my $value (keys $max_set_output) {
      $header{$value} = $max_set_output->{$value};
    }
  }
  elsif($self->{config_data}) {
    foreach my $value (keys $self->{config_data_hash}) {
      $header{$value} = $self->{config_data_hash}->{$value};
    }
  }
  else {
    foreach my $value (keys %$default_output) {
      $header{$value} = $default_output->{$value};
    }
  }

  return \%header;
}

sub run {
  my ($self, $tva) = @_;

  return {} unless grep {$_->SO_term =~ 'missense_variant'} @{$tva->get_all_OverlapConsequences};
  my $tv = $tva->transcript_variation;

  my $peptide_start = defined($tv->translation_start) ? $tv->translation_start : undef;
  my $peptide_end = defined($tv->translation_end) ? $tv->translation_end : undef;
  my $transcript_id = $tva->transcript->translation->stable_id;

  return {} unless defined($transcript_id) && defined($peptide_start) && defined($peptide_end);

  my $query = $self->get_sparql_query($peptide_start,$transcript_id);

  # run SPARQL query
  my $query_output;
  eval {
    $query_output = `curl -X POST -H "Accept:application/sparql-results+json" --data-urlencode "query=$query" https://sparql.nextprot.org/ 2> /dev/null`;
  };
  warn $@ if $@;

  my $output = decode_json ($query_output);

  my %result_hash;
  my %result_hash_final;

  # Output format: 'iso','spos','epos','annot_type','callret-4'
  # 'iso' -> isoform URL to neXtProt page; 'spos' -> start position; 'epos' -> end position; 'annot_type' -> annotation type (e.g. PdbMapping, Variant, etc.);
  # 'callret-4' -> data
  my $output_list = $output->{results}->{bindings};
  return {} if (@$output_list == 0);

  foreach my $results (@$output_list) {
    my $isoform_url = $results->{iso}->{value};
    my $start_pos = $results->{spos}->{value};
    my $end_pos = $results->{epos}->{value};
    my $annot_type = $results->{annot_type}->{value};
    $annot_type =~ s/.*#//;
    my $data = $results->{'callret-4'}->{value};
    # PdbMapping and Variant values contain ';'
    if($data =~ /;/) {
      if($annot_type eq 'Variant') { 
        $data =~ s/;/\./g; 
      }
      else {
        $data =~ s/;//g;
      }
    }
    $data =~ s/\.$//;

    # There is only one URL
    if($self->{url} && !$result_hash{'neXtProt_url'}) {
      my @isoform_value = ($isoform_url);
      $result_hash{'neXtProt_url'} = \@isoform_value;
    }

    # Some annot_type have more than one value
    # Need to check if it's not duplicated
    if($result_hash{'neXtProt_'.$annot_type}) {
      my $annot_type_data = $start_pos.','.$end_pos.','.$data;
      push @{$result_hash{'neXtProt_'.$annot_type}}, $annot_type_data unless grep{$_ eq $annot_type_data} @{$result_hash{'neXtProt_'.$annot_type}};
    }
    else {
      my @list_of_data;
      push @list_of_data, $start_pos.','.$end_pos.','.$data;
      $result_hash{'neXtProt_'.$annot_type} = \@list_of_data;
    }
  }

  my @keys;
  if($self->{max_set}) {
    @keys = keys %$max_set_output;
  }
  elsif($self->{config_data}) {
    @keys = keys $self->{config_data_hash};
  }
  else {
    @keys = keys %$default_output;
  }

  foreach my $key (@keys) {
    if($result_hash{$key}) {
      my $data_to_return = $result_hash{$key};
      my $join_data = join('|', @$data_to_return);
      $result_hash_final{$key} = $join_data;
    }
    elsif(!$result_hash{$key} && $self->{all_labels}) {
      $result_hash_final{$key} = '-';
    }
  }

  return \%result_hash_final;
}

sub get_sparql_query {
  my ($self, $peptide_start, $transcript_id) = @_;

  my $query = "PREFIX : <http://nextprot.org/rdf#>
               PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
               PREFIX up: <http://purl.uniprot.org/core/>
               PREFIX isoform: <http://nextprot.org/rdf/isoform/>
               select distinct ?iso ?spos ?epos ?annot_type str(?txt)
               where {
                 values ?poi {$peptide_start}
                 values ?ensp {'$transcript_id'} # for e101 mapping is for canoncila - somtimes our cnonincal is different than the uniprot canonical
                 bind (IRI(CONCAT('http://rdf.ebi.ac.uk/resource/ensembl.protein/',?ensp)) as ?ENSP_IRI)
                 SERVICE <http://sparql.uniprot.org/sparql> {
                   SELECT * WHERE {
                     ?enst up:translatedTo ?ENSP_IRI .
                     ?enst rdfs:seeAlso  ?upiso .
                   }
                 }
               BIND(IRI(replace(str(?upiso),'http://purl.uniprot.org/isoforms/','http://nextprot.org/rdf/isoform/NX_')) AS ?iso) .
                 ?entry :isoform ?iso .
                 ?iso :positionalAnnotation ?statement .
                 ?statement rdfs:comment ?txt .
                 ?statement a ?annot_type .
                 ?statement :start ?spos; :end ?epos .
                 filter((?spos <= ?poi) && (?epos >= ?poi))
               } order by ?spos";

  return $query;
}

sub build_data_hash {
  my $self = shift;
  my $file = 'neXtProt_headers.txt';

  my %headers_file_hash;
  my %output_hash;

  if (! -e $file) {
    die ("ERROR: neXtProt_headers file is not available.");
  } else {
    open FILE, $file;
    while(<FILE>) {
      chomp;
      my ($value, $description) = split(/\t/, $_);
      die ("ERROR: neXtProt value is missing from file.") if(!$value);
      $headers_file_hash{$value} = $description;
    }
    close FILE;
  }

  my @data_from_user = split(/[\;\&\|]/, $self->{config_data});
  foreach my $dfu (@data_from_user) {
    if($headers_file_hash{$dfu}) {
      $self->{config_data_hash}->{'neXtProt_' . $dfu} = $headers_file_hash{$dfu};
    }
    else {
      die ("ERROR: $dfu is not available in neXtProt. Check file 'neXtProt_headers.txt' to see the data that is valid to query.\n");
    }
  }
}

1;

