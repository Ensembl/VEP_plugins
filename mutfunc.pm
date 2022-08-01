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

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);


my @ALL_AAS = qw(A C D E F G H I K L M N P Q R S T V W Y);

# dispatch table for frunctions that parse item values
my $parse_item_value = {
  motif => \&parse_motif,
  tfbs  => \&parse_tfbs
};

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  
  my $param_hash = $self->params_to_hash();

  die "ERROR: db is not specified which is a mandatory parameter\n" unless defined $param_hash->{db};

  $self->{db} = $param_hash->{db};

  $self->{motif} = 1 if $param_hash->{motif} || $param_hash->{all};
  $self->{int} = 1 if $param_hash->{int} || $param_hash->{all};
  $self->{mod} = 1 if $param_hash->{mod} || $param_hash->{all};
  $self->{exp} = 1 if $param_hash->{exp} || $param_hash->{all};

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
  "lost - '1' if the mutation causes the motif to be lost and '0' otherwise" if defined $self->{motif};

  $header{"mutfunc_int"} = "Interaction interfaces destabilization analysis from mutfunc db. ".
  "Output fields are seperated by comma and include: ".
  "evidence - 'EXP' for experimental model and 'MDL' for homology models and 'MDD' for domain-domain homology models, ".
  "dG_wt - reference interface energy (kcal/mol), ".
  "dG_mt - mutated interface energy (kcal/mol), ".
  "ddG - change in interface stability between mutated and reference structure (kcal/mol) mutations where ddG >= 2 kcal/mol can be considered deleterious, ".
  "dG_wt_sd - dG_wt standard deviation (kcal/mol), ".
  "dG_mt_sd - dG_mt standard deviation (kcal/mol), ".
  "ddG_sd - ddG standard deviation (kcal/mol), " if defined $self->{int};

  $header{"mutfunc_mod"} = "Protein structure destabilization analysis (homology models) from mutfunc db. ".
  "Output fields are seperated by comma and include: ".
  "dG_wt - reference structure energy (kcal/mol), ".
  "dG_mt - mutated structure energy (kcal/mol), ".
  "ddG - change in structure stability between mutated and reference structure (kcal/mol) mutations where ddG >= 2 kcal/mol can be considered deleterious, ".
  "dG_wt_sd - dG_wt standard deviation (kcal/mol), ".
  "dG_mt_sd - dG_mt standard deviation (kcal/mol), ".
  "ddG_sd - ddG standard deviation (kcal/mol), " if defined $self->{mod};

  $header{"mutfunc_exp"} = "Protein structure destabilization analysis (experimental models) from mutfunc db. ".
  "Output fields are seperated by comma and include: ".
  "dG_wt - reference structure energy (kcal/mol), ".
  "dG_mt - mutated structure energy (kcal/mol), ".
  "ddG - change in structure stability between mutated and reference structure (kcal/mol) mutations where ddG >= 2 kcal/mol can be considered deleterious, ".
  "dG_wt_sd - dG_wt standard deviation (kcal/mol), ".
  "dG_mt_sd - dG_mt standard deviation (kcal/mol), ".
  "ddG_sd - ddG standard deviation (kcal/mol), " if defined $self->{exp};
   
  return \%header;
}

sub expand_matrix {
  my ($matrix) = @_;
  my $expanded_matrix = Compress::Zlib::memGunzip($matrix) or 
    throw("Failed to gunzip: $gzerrno");

  return $expanded_matrix;
}

sub retrieve_item_value {
  my ($matrix, $pos, $aa, $tot_packed_len) = @_;

  my $item_value = (substr $matrix, $pos * $tot_packed_len * 20 + $aa * $tot_packed_len, $tot_packed_len );

  return $item_value;
}

sub parse_motif {
  my ($item_value) = @_;

  my ($elm, $lost) = unpack "A24v", $item_value;

  $elm = undef if $elm eq "undefined";
  $lost = undef if $lost == 0xFFFF;

  return $elm, $lost;
}

sub parse_destabilizers {
  my ($item_value, $item) = @_;
  my @evidence_lval = qw(EXP MDD MDL);
  my $evidence_val;

  # only int item type have evidence
  if ($item eq "int"){
    my $evidence = unpack "v", $item_value;

    # now omit the evidence part from item value
    $item_value = substr $item_value, 2; 

    # get the value of the evidence
    $evidence = undef $evidence == 0xFFFF;
    my $evidence_val = defined $evidence ? $evidence_lval[$evidence] : undef;
  }
  
  # get the rest
  my ($dG_wt, $ddG, $dG_wt_sd, $dG_mt_sd, $ddG_sd) = unpack "A8A8A8A8A8", $item_value;

  $dG_wt = undef if $dG_wt eq "10000000";
  $ddG = undef if $ddG eq "10000000";
  $dG_wt_sd = undef if $dG_wt_sd eq "10000000";
  $dG_mt_sd = undef if $dG_mt_sd eq "10000000";
  $ddG_sd = undef if $ddG_sd eq "10000000";

  return $evidence_val, $dG_wt, $ddG, $dG_wt_sd, $dG_mt_sd, $ddG_sd if $item eq "int";
  return $dG_wt, $ddG, $dG_wt_sd, $dG_mt_sd, $ddG_sd;
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

    # get matrix for each item value and parse them
    if ($item eq "motif" && $self->{motif}) {
      my $tot_packed_len = 26;
      my $item_value = retrieve_item_value($expanded_matrix, $pos, $peptide_number, $tot_packed_len);

      if ($item_value){
        my ($elm, $lost) = $parse_item_value->{$item}($item_value);

        $result->{mutfunc_motif} = $elm if defined $elm;
        $result->{mutfunc_motif} .= defined $lost ? ",".$lost : ",";

        # delete if there is no defined value
        delete $result->{mutfunc_motif} unless $result->{mutfunc_motif} =~ /[^,]/;
      }
    }
    elsif ($item eq "int" && $self->{int}) {
      my $tot_packed_len = 42;
      my $item_value = retrieve_item_value($expanded_matrix, $pos, $peptide_number, $tot_packed_len);

      if ($item_value){
        my ($evidence, $dG_wt, $ddG, $dG_wt_sd, $dG_mt_sd, $ddG_sd) = parse_destabilizers($item_value, $item);

        # dG_mt can be calculated from dG_wt and ddG
        my $dG_mt = (defined $dG_wt && defined $ddG) ? $dG_wt + $ddG : undef;

        $result->{mutfunc_int} = $evidence if defined $evidence;
        $result->{mutfunc_int} .= defined $dG_wt ? ",".$dG_wt : ",";
        $result->{mutfunc_int} .= defined $dG_mt ? ",".$dG_mt : ",";
        $result->{mutfunc_int} .= defined $ddG ? ",".$ddG : ",";
        $result->{mutfunc_int} .= defined $dG_wt_sd ? ",".$dG_wt_sd : ",";
        $result->{mutfunc_int} .= defined $dG_mt_sd ? ",".$dG_mt_sd : ",";
        $result->{mutfunc_int} .= defined $ddG_sd ? ",".$ddG_sd : ",";

        # delete if there is no defined value
        delete $result->{mutfunc_int} unless $result->{mutfunc_int} =~ /[^,]/;
      }
    }
    elsif ($item eq "mod" || $item eq "exp") {
      if ($self->{$item}) {
        my $tot_packed_len = 40;
        my $item_value = retrieve_item_value($expanded_matrix, $pos, $peptide_number, $tot_packed_len);

        if ($item_value){
          my ($dG_wt, $ddG, $dG_wt_sd, $dG_mt_sd, $ddG_sd) = parse_destabilizers($item_value, $item);;

          # dG_mt can be calculated from dG_wt and ddG
          my $dG_mt = (defined $dG_wt && defined $ddG) ? $dG_wt + $ddG : undef;

          my $hash_key = "mutfunc_" . $item;
          $result->{$hash_key} = $dG_wt if defined $dG_wt;
          $result->{$hash_key} .= defined $dG_mt ? ",".$dG_mt : ",";
          $result->{$hash_key} .= defined $ddG ? ",".$ddG : ",";
          $result->{$hash_key} .= defined $dG_wt_sd ? ",".$dG_wt_sd : ",";
          $result->{$hash_key} .= defined $dG_mt_sd ? ",".$dG_mt_sd : ",";
          $result->{$hash_key} .= defined $ddG_sd ? ",".$ddG_sd : ",";

          # delete if there is no defined value
          delete $result->{$hash_key} unless $result->{$hash_key} =~ /[^,]/;
        }
      }
    }
    
  }

  return $result;
}

1;