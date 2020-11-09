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
 retrieves data for missense and stop gain variants from neXtProt, which is a comprehensive 
 human-centric discovery platform that offers integration of and navigation 
 through protein-related data for example, variant information, localization 
 and interactions  (https://www.nextprot.org/).

 Please cite the neXtProt publication alongside the VEP if you use this resource:
 https://doi.org/10.1093/nar/gkz995

 This plugin is only suitable for small sets of variants as an additional 
 individual remote API query is run for each variant.

 Running options:
 (Default) the data retrieved by default is the MatureProtein, NucleotidePhosphateBindingRegion,
 Variant, MiscellaneousRegion, TopologicalDomain and InteractingRegion.
 The plugin can also be run with other options to retrieve other data than the default.
 
 Options are passed to the plugin as key=value pairs:
 max_set        : Set value to 1 to return all available protein-related data 
                  (includes the default data)

 return_values  : The set of data to be returned with different data separated by '&'. 
                  Use file 'neXtProt_headers.txt' to check which data (labels) are available.
                  Example: --plugin neXtProt,return_values='Domain&InteractingRegion'

 url            : Set value to 1 to include the URL to link to the neXtProt entry.

 all_labels     : Set value to 1 to include all labels, even if data is not available.

 position       : Set value to 1 to include the start and end position in the protein.

 * note: 'max_set' and 'return_values' cannot be used simultaneously.


 Output:
  By default, the plugin only returns data that is available. Example (default behaviour):
  neXtProt_MatureProtein=Rho guanine nucleotide exchange factor 10

  The option 'all_labels' returns a consistent set of the requested fields, using "-" where 
  values are not available. Same example as above:
  neXtProt_MatureProtein=Rho guanine nucleotide exchange factor 10;
  neXtProt_InteractingRegion=-;neXtProt_NucleotidePhosphateBindingRegion=-;neXtProt_Variant=-;
  neXtProt_MiscellaneousRegion=-;neXtProt_TopologicalDomain=-;

  Of notice, multiple values can be returned for the same label. In this case, the values will
  be separeted by '|' for tab and txt format, and '&' for VCF format. 

 The plugin can then be run as default:
 ./vep -i variations.vcf --plugin neXtProt

 or to return only the data specified by the user:
 ./vep -i variations.vcf --plugin neXtProt,return_values='Domain&InteractingRegion'


=cut

package neXtProt;

use strict;
use warnings;
use JSON::XS;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my $default_output = {
  'neXtProt_MatureProtein'                    => 'Extent of an active peptide or a polypetide chain in the mature protein',
  'neXtProt_NucleotidePhosphateBindingRegion' => 'Nucleotide phosphate binding region',
  'neXtProt_Variant'                          => 'Variant-specific annotations',
  'neXtProt_MiscellaneousRegion'              => 'Region of interest in the sequence',
  'neXtProt_TopologicalDomain'                => 'Location of non-membrane regions of membrane-spanning proteins',
  'neXtProt_InteractingRegion'                => 'Region interacting with another macromolecule'
};

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $param_hash = $self->params_to_hash();

  if(defined($param_hash->{max_set}) && defined($param_hash->{return_values})) {
    die "ERROR: Can't use max_set and return_values simultaneously!\n";
  }

  # Return the isoform URL to the neXtProt web page
  if(defined($param_hash->{url})) {
    $self->{url} = $param_hash->{url};

    if(defined($param_hash->{max_set}) || defined($param_hash->{return_values})) {
      $self->{return_values_hash}->{'neXtProt_url'} = 'neXtProt URL';
    }
    else {
      $default_output->{'neXtProt_url'} = 'neXtProt URL';
    }
  }

  if(defined($param_hash->{max_set})) {
    $self->{max_set} = $param_hash->{max_set};
    $self->build_data_hash(1);
  }

  if(defined($param_hash->{return_values})) {
    $self->{return_values} = $param_hash->{return_values};
    $self->build_data_hash(0);
  }

  if(defined($param_hash->{all_labels})) {
    $self->{all_labels} = $param_hash->{all_labels};
  }

  if(defined($param_hash->{position})) {
    $self->{position} = $param_hash->{position};
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
    foreach my $value (keys %{$self->{return_values_hash}}) {
      $header{$value} = $self->{return_values_hash}->{$value};
    }
  }
  elsif($self->{return_values}) {
    foreach my $value (keys %{$self->{return_values_hash}}) {
      $header{$value} = $self->{return_values_hash}->{$value};
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

  return {} unless grep {$_->SO_term =~ 'missense_variant|stop_gain'} @{$tva->get_all_OverlapConsequences};
  my $tv = $tva->transcript_variation;

  my $peptide_start = defined($tv->translation_start) ? $tv->translation_start : undef;
  my $translation_id = $tva->transcript->translation->stable_id;

  my $tl = $tva->transcript->translation;
  my @uniprot = @{$tl->get_all_DBLinks('Uniprot_isoform')};

  return {} unless defined($peptide_start) && scalar @uniprot > 0; 
  
  my $isoform_id = $uniprot[0]->display_id;

  my $query = $self->get_sparql_query($peptide_start,$isoform_id);

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
  # 'nx_lnk' -> isoform URL to neXtProt page; 'spos' -> start position; 'epos' -> end position; 'annot_type' -> annotation type (e.g. PdbMapping, Variant, etc.);
  # 'annot_descr' -> data
  my $output_list = $output->{results}->{bindings};
  return {} if (@$output_list == 0);

  foreach my $results (@$output_list) {

    my $isoform_url = $results->{nx_lnk}->{value};
    my $start_pos = $results->{spos}->{value};
    my $end_pos = $results->{epos}->{value};
    my $annot_type = $results->{annot_type}->{value};
    $annot_type =~ s/.*#//;
    my $data = $results->{annot_descr}->{value};
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
      my $annot_type_data = $self->{position} ? $start_pos.','.$end_pos.','.$data : $data;
      push @{$result_hash{'neXtProt_'.$annot_type}}, $annot_type_data unless grep{$_ eq $annot_type_data} @{$result_hash{'neXtProt_'.$annot_type}};
    }
    else {
      my @list_of_data;
      my $annot_type_data = $self->{position} ? $start_pos.','.$end_pos.','.$data : $data;
      push @list_of_data, $annot_type_data;
      $result_hash{'neXtProt_'.$annot_type} = \@list_of_data;
    }
  }

  my @keys;
  if($self->{max_set} || $self->{return_values}) {
    @keys = keys %{$self->{return_values_hash}};
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
  my ($self, $peptide_start, $translation_id) = @_;

  my $query = "PREFIX : <http://nextprot.org/rdf#>
               PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
               PREFIX up: <http://purl.uniprot.org/core/>
               PREFIX isoform: <http://nextprot.org/rdf/isoform/>
               PREFIX proteoform: <http://nextprot.org/rdf/proteoform/>
               PREFIX cv: <http://nextprot.org/rdf/terminology/>
               select distinct ?iso ?nx_lnk ?spos ?epos ?annot_type str(?txt) as ?annot_descr where {
                 values ?upiso {'$translation_id'}
                 BIND(IRI(CONCAT('http://nextprot.org/rdf/isoform/NX_',?upiso)) AS ?iso) .
                 ?entry :isoform ?iso .
                 values ?poi {$peptide_start}
                 {
                   ?iso :positionalAnnotation ?statement .
                   ?statement rdfs:comment ?txt .
	           bind(replace(str(?iso),'http://nextprot.org/rdf/isoform/','') as  ?iso_ac) .
                   bind(concat('https://www.nextprot.org/entry/', ?iso_ac, '/sequence?isoform=', ?iso_ac) as ?nx_lnk) .
                   ?statement a ?annot_type .
                   ?statement :start ?spos; :end ?epos .
                 }
                 union
                 {
                   ?iso :proteoform ?pf .
                   ?pf :modification ?varmut; :phenotypicVariation ?phvar .
                   ?varmut :start ?spos; :end ?epos.
                   ?phvar :term ?phtype; :impactedObject / :term / rdfs:label ?ioTermlab .
                   ?phvar a ?annot_type; :entryAnnotationId ?eid .
                   ?phtype :childOf cv:ME_0000002; rdfs:label ?effect .
                   bind (concat(CONCAT(?effect,' '),?ioTermlab) as ?txt)
                 }
                 filter((?spos <= ?poi) && (?epos >= ?poi))
                 } order by ?spos";

  return $query;
}

sub build_data_hash {
  my $self = shift;
  my $option = shift; # Set to 1 to return all data from header file, set to 0 to return data specified by the user

  my $plugin_dir = $INC{'neXtProt.pm'};
  $plugin_dir =~ s/neXtProt\.pm//i;
  my $file = $plugin_dir.'/neXtProt_headers.txt';

  my %headers_file_hash;
  my %output_hash;

  if (! -e $file) {
    die ("ERROR: neXtProt_headers file is not available in $plugin_dir");
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

  # Going to return the data specified by the user
  if($option == 0) {
    my @data_from_user = split(/[\;\&\|]/, $self->{return_values});
    foreach my $dfu (@data_from_user) {
      if($headers_file_hash{$dfu}) {
        $self->{return_values_hash}->{'neXtProt_' . $dfu} = $headers_file_hash{$dfu};
      }
      else {
        die ("ERROR: $dfu is not available in neXtProt. Check file 'neXtProt_headers.txt' to see the data that is valid to query.\n");
      }
    }
  }
  # Return all data available in header file
  else {
    foreach my $hf (keys %headers_file_hash) {
      $self->{return_values_hash}->{'neXtProt_' . $hf} = $headers_file_hash{$hf};
    }
  }
}

1;

