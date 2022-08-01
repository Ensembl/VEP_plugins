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

 mutfunc

=head1 SYNOPSIS

 mv mutfunc.pm ~/.vep/Plugins
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
 
 By default the output will always contain feature_type and interaction_ac from the IntAct data file. You can also add more fields using the following options -
 feature_ac			: Set value to 1 to include Feature AC in the output
 feature_short_label		: Set value to 1 to include Feature short label in the output
 feature_annotation		: Set value to 1 to include Feature annotation in the output
 ap_ac				: Set value to 1 to include Affected protein AC in the output
 interaction_participants	: Set value to 1 to include Interaction participants in the output
 pmid				: Set value to 1 to include PubMedID in the output

 There are also two other options for customizing the output - 
 all                            : Set value to 1 to include all the fields
 minimal                        : Set value to 1 to overwrite default behavior and include only interaction_ac 
				  in the output by default

 See what this options mean - https://www.ebi.ac.uk/intact/download/datasets#mutations
 
 Note that, interaction accession can be used to link to full details on the interaction website. For example, 
 where the VEP output reports an interaction_ac of EBI-12501485, the URL would be : 

   https://www.ebi.ac.uk/intact/details/interaction/EBI-12501485

=cut

package mutfunc;

use strict;
use warnings;
use DBI;
use Compress::Zlib;
use Digest::MD5 qw(md5_hex);
use List::MoreUtils qw(first_index);

use Bio::EnsEMBL::Utils::Exception qw(warning);
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);


my @ALL_AAS = qw(A C D E F G H I K L M N P Q R S T V W Y);

# dispatch table for frunctions to retrieve items
my $retrieve_item = {
  motif_elm => \&retrieve_motif_elm,
  motif_lost => \&retrieve_motif_lost,
  int_evidence => \&retrieve_int_evidence,
  energy => \&retrieve_energy
};

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  
  my $param_hash = $self->params_to_hash();

  die "ERROR: db is not specified which is a mandatory parameter\n" unless defined $param_hash->{db};

  $self->{db} = $param_hash->{db};

  $self->{motif} = 1 if $param_hash->{motif} || $param_hash->{all};
  $self->{int} = 1 if $param_hash->{int} || $param_hash->{all};

  $self->{initial_pid} = $$;

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;
  
  my %header;

  $header{"mutfunc_motif"} = "Nonsynonymous mutations impact on linear motif from mutfunc db. ".
  "Output fields are seperated by comma and include: ".
  "elm - ELM accession of the linear motif, ".
  "lost - '1' if the mutation causes the motif to be lost and '0' otherwise";

  $header{"mutfunc_int"} = "Interaction interfaces destabilization analysis from mutfunc db. ".
  "Output fields are seperated by comma and include: ".
  "evidence - 'EXP' for experimental model and 'MDL' for homology models and 'MDD' for domain-domain homology models, ".
  "dG_wt - reference interface energy (kcal/mol), ".
  "dG_mt - mutated interface energy (kcal/mol), ".
  "ddG - change in interface stability between mutated and reference structure (kcal/mol) mutations with ddG >= 2 kcal/mol can be considered deleterious, ".
  "dG_wt_sd - dG_wt standard deviation (kcal/mol), ".
  "dG_mt_sd - dG_mt standard deviation (kcal/mol), ".
  "ddG_sd - ddG standard deviation (kcal/mol), ";
   
  return \%header;
}

sub expand_matrix {
  my ($matrix) = @_;

  my $expanded_matrix = Compress::Zlib::memGunzip($matrix) or 
    throw("Failed to gunzip: $gzerrno");

  return $expanded_matrix;
}

sub retrieve_motif_elm {
  my ($matrix, $pos, $aa) = @_;

  my $packed_len = 24;
  my $item_value = (unpack "A24", substr($matrix, $pos * $packed_len * 20 + $aa * $packed_len, $packed_len) );

  return undef if $item_value eq "undefined";
  return $item_value;
}

sub retrieve_motif_lost {
  my ($matrix, $pos, $aa) = @_;

  my $packed_len = 2;
  my $item_value = (unpack "v", substr($matrix, $pos * $packed_len * 20 + $aa * $packed_len, $packed_len) );

  return undef if $item_value == 0xFFFF;
  return $item_value;
}

sub retrieve_int_evidence {
  my ($matrix, $pos, $aa) = @_;

  my @evidence_sval = qw(EXP MDD MDL);

  my $packed_len = 2;
  my $item_value = (unpack "v", substr($matrix, $pos * $packed_len * 20 + $aa * $packed_len, $packed_len) );

  return undef if $item_value == 0xFFFF;
  return $evidence_sval[$item_value];
}

sub retrieve_energy {
  my ($matrix, $pos, $aa) = @_;

  my $packed_len = 8;
  my $item_value = (unpack "A8", substr($matrix, $pos * $packed_len * 20 + $aa * $packed_len, $packed_len) );

  return undef if $item_value eq "100000000";
  return $item_value;
}

sub run {
  my ($self, $tva) = @_;

  # get the trascript related to the variant
  my $tr = $tva->transcript;

  # get the translation
  my $translation = $tr->translate;
  return {} unless $translation;

  # get the md5 hash of the peptide sequence
  my $md5 = md5_hex($translation->seq);

  # forked, reconnect to DB
  if($$ != $self->{initial_pid}) {
    $self->{dbh} = DBI->connect("dbi:SQLite:dbname=".$self->{db},"","");
    $self->{get_sth} = $self->{dbh}->prepare("SELECT md5, item, matrix FROM consequences WHERE md5 = ?");

    # set this so only do once per fork
    $self->{initial_pid} = $$;
  }
  $self->{get_sth}->execute($md5);

  my $result = {};
  while (my $arrayref = $self->{get_sth}->fetchrow_arrayref) {
    my $item = $arrayref->[1];
    my $matrix = $arrayref->[2];

    # expand the compressed matrix
    my $expanded_matrix = expand_matrix($matrix);

    # position and peptide to retrieve value from matrix
    my $pos = $tva->transcript_variation->translation_start;
    my $peptide = $tva->peptide;

    # in matrix position is 0 indexed
    $pos--;

    # we need position of peptide in the ALL_AAS array
    my $peptide_number = (first_index { $_ eq $peptide } @ALL_AAS);

    if ($item =~ /motif*/ && $self->{motif}) {
      my $item_value = $retrieve_item->{$item}($expanded_matrix, $pos, $peptide_number);

      if ($item_value){
        unless (defined $result->{mutfunc_motif}) {
          $result->{mutfunc_motif} = $item_value;
        }
        else{
          $result->{mutfunc_motif} = $result->{mutfunc_motif}.",".$item_value;
        }
      }
    }
    elsif ($item =~ /int*/ && $self->{int}) {
      my $item_retrieval = $item eq "int_evidence" ? $item : "energy";
      my $item_value = $retrieve_item->{$item_retrieval}($expanded_matrix, $pos, $peptide_number);

      if ($item_value){
        unless (defined $result->{mutfunc_int}) {
          $result->{mutfunc_int} = $item_value;
        }
        else{
          $result->{mutfunc_int} = $result->{mutfunc_int}.",".$item_value;
        }
      }
    }
  }

  return $result;
}

1;