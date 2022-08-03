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

use Data::Dumper;
use strict;
use warnings;
use DBI;
use Compress::Zlib;
use Digest::MD5 qw(md5_hex);
use List::MoreUtils qw(first_index);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my @ALL_AAS = qw(A C D E F G H I K L M N P Q R S T V W Y);

my $field_order = {
  motif => ["elm", "lost"],
  int   => ["evidence", "dG_wt", "ddG", "dG_wt_sd", "dG_mt_sd", "ddG_sd"],
  mod   => ["dG_wt", "ddG", "dG_wt_sd", "dG_mt_sd", "ddG_sd"],
  exp   => ["dG_wt", "ddG", "dG_wt_sd", "dG_mt_sd", "ddG_sd"],
  tfbs  => ["impact", "tf", "downstream", "g_strand", "s_strand", "wt_score", "mt_score", "ic_diff", "cells"]
};

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  
  my $param_hash = $self->params_to_hash();

  $self->{tfbs} = 1 if $param_hash->{tfbs} || $param_hash->{all};

  die "ERROR: tfbs file not specified but tfbs output is enabled\n" if ( (defined $self->{tfbs}) && !(defined $param_hash->{file}) );
  $self->{file} = $param_hash->{file};

  $self->add_file($self->{file}) if defined $self->{file};

  $self->{motif} = 1 if $param_hash->{motif} || $param_hash->{all};
  $self->{int} = 1 if $param_hash->{int} || $param_hash->{all};
  $self->{mod} = 1 if $param_hash->{mod} || $param_hash->{all};
  $self->{exp} = 1 if $param_hash->{exp} || $param_hash->{all};

  die "ERROR: db is not specified but some of the options enabled require it\n" if ( 
    ( (defined $self->{motif}) || 
      (defined $self->{int}) || 
      (defined $self->{mod}) || 
      (defined $self->{exp}) 
    ) && 
    !(defined $param_hash->{db}) 
  );
  $self->{db} = $param_hash->{db};

  if( ($self->{config}->{output_format} eq "json") || $self->{config}->{rest}){
    $self->{output_json} = 1;
  }

  $self->{initial_pid} = $$;

  return $self;
}

sub feature_types {
  return ['Transcript', 'Intergenic'];
}

sub get_header_info {
  my ($self) = shift;
  
  my %header;

  $header{"mutfunc_motif"} = "Nonsynonymous mutations impact on linear motif from mutfunc db. ".
  "Output fields are separated by ',' (or '&' for vcf format) and include: ".
  "elm - ELM accession of the linear motif, ".
  "lost - '1' if the mutation causes the motif to be lost and '0' otherwise" if defined $self->{motif};

  $header{"mutfunc_int"} = "Interaction interfaces destabilization analysis from mutfunc db. ".
  "Output fields are separated by ',' (or '&' for vcf format) and include: ".
  "evidence - 'EXP' for experimental model and 'MDL' for homology models and 'MDD' for domain-domain homology models, ".
  "dG_wt - reference interface energy (kcal/mol), ".
  "dG_mt - mutated interface energy (kcal/mol), ".
  "ddG - change in interface stability between mutated and reference structure (kcal/mol) mutations where ddG >= 2 kcal/mol can be considered deleterious, ".
  "dG_wt_sd - dG_wt standard deviation (kcal/mol), ".
  "dG_mt_sd - dG_mt standard deviation (kcal/mol), ".
  "ddG_sd - ddG standard deviation (kcal/mol), " if defined $self->{int};

  $header{"mutfunc_mod"} = "Protein structure destabilization analysis (homology models) from mutfunc db. ".
  "Output fields are separated by ',' (or '&' for vcf format) and include: ".
  "dG_wt - reference structure energy (kcal/mol), ".
  "dG_mt - mutated structure energy (kcal/mol), ".
  "ddG - change in structure stability between mutated and reference structure (kcal/mol) mutations where ddG >= 2 kcal/mol can be considered deleterious, ".
  "dG_wt_sd - dG_wt standard deviation (kcal/mol), ".
  "dG_mt_sd - dG_mt standard deviation (kcal/mol), ".
  "ddG_sd - ddG standard deviation (kcal/mol), " if defined $self->{mod};

  $header{"mutfunc_exp"} = "Protein structure destabilization analysis (experimental models) from mutfunc db. ".
  "Output fields are separated by ',' (or '&' for vcf format) and include: ".
  "dG_wt - reference structure energy (kcal/mol), ".
  "dG_mt - mutated structure energy (kcal/mol), ".
  "ddG - change in structure stability between mutated and reference structure (kcal/mol) mutations where ddG >= 2 kcal/mol can be considered deleterious, ".
  "dG_wt_sd - dG_wt standard deviation (kcal/mol), ".
  "dG_mt_sd - dG_mt standard deviation (kcal/mol), ".
  "ddG_sd - ddG standard deviation (kcal/mol), " if defined $self->{exp};

  $header{"mutfunc_tfbs"} = "Transcription binding sites disruption analysis from mutfunc db. ".
  "Output fields are separated by ',' (or '&' for vcf format) and include: ".
  "impact -  1 if mutation is expected to disrupt the TFBS and 0 otherwise, ".
  "tf - transcription factor predicted to bind this binding site, ".
  "downstream - gene downstream from this TFBS, ".
  "g_strand - strand containing the downstream gene, ".
  "s_strand - strand containing binding site, ".
  "wt_score - binding score of the transcription factor to the wildtype site, ".
  "mt_score - binding score of the transcription factor to the mutated site".
  "ic_diff - information content difference between the wildtype and mutant bases".
  "cells - cell lines or tissues with evidence of chip-seq (human only) multiple cells are separated by 'and'" if defined $self->{tfbs};
   
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
    $evidence = undef if $evidence == 0xFFFF;
    $evidence_val = defined $evidence ? $evidence_lval[$evidence] : undef;
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

sub format_output{
  my ($self, $data, $item) = @_;

  my $result = {};

  if ($self->{output_json}){
    my %hash = map { $_ => $data->{$_} } @{ $field_order->{$item} };
    $result->{$item} = \%hash;
  }
  else{
    my $key = "mutfunc_" . $item;
    $result->{$key} = join(",", map { $data->{$_} } @{ $field_order->{$item} });
  }

  return $result;
}

sub process_from_file {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  
  # get allele
  my $allele = $tva->variation_feature_seq;
  return {} unless $allele =~ /^[ACGT-]+$/;

  # get matched line from the file
  my @data =  @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})};

  # parse data and generate output
  my $result_from_file = {};
  foreach (@data) {
    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => [$allele],
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
       ref  => $_->{ref},
       alts => [$_->{alt}],
       pos  => $_->{start},
      }
    );

    if (@$matches){
      my $result = $_->{result}; 
      
      $result->{cells} =~ s/,/and/g;

      # we are not expecting multiple results so no code for joiing multiple results
      $result_from_file = $self->format_output($result, "tfbs");
    }
  }

  return $result_from_file;
}

sub process_from_db {
  my ($self, $tva) = @_;

  # get the trascript related to the variant
  return {} unless $tva->can('transcript');
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

  my $result_from_db = {};
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
    # motif
    if ($item eq "motif" && $self->{motif}) {
      # get the item value for specific pos and amino acid from matrix
      my $tot_packed_len = 26;
      my $item_value = retrieve_item_value($expanded_matrix, $pos, $peptide_number, $tot_packed_len);

      if ($item_value){
        # parse the item value from matrix
        my ($elm, $lost) = parse_motif($item_value);

        # format the output
        my $formatted_output = {};
        if(defined $elm || defined $lost){
          $formatted_output = $self->format_output({
            elm   => $elm,
            lost  => $lost
          }, $item);
        }

        @$result_from_db{ keys %$formatted_output } = values %$formatted_output;
      }
    }
    # int and mod and exp
    elsif ($item eq "int" || $item eq "mod" || $item eq "exp") {
      if ($self->{$item}){
        # get the item value for specific pos and amino acid from matrix
        my $tot_packed_len = ($item eq "int") ? 42 : 40;
        my $item_value = retrieve_item_value($expanded_matrix, $pos, $peptide_number, $tot_packed_len);

        if ($item_value){
          # parse the item value from matrix
          my ($evidence, $dG_wt, $ddG, $dG_wt_sd, $dG_mt_sd, $ddG_sd);
          if ($item eq "int"){
            ($evidence, $dG_wt, $ddG, $dG_wt_sd, $dG_mt_sd, $ddG_sd) = parse_destabilizers($item_value, $item);
          }
          else{
            ($dG_wt, $ddG, $dG_wt_sd, $dG_mt_sd, $ddG_sd) = parse_destabilizers($item_value, $item);
          }

          # dG_mt can be calculated from dG_wt and ddG
          my $dG_mt = (defined $dG_wt && defined $ddG) ? $dG_wt + $ddG : undef;

          # format the output
          my $formatted_output = {};
          if(defined $evidence || defined $dG_wt || defined $ddG || defined $dG_wt_sd || defined $dG_mt_sd || defined $ddG_sd){
            my $data = {
              dG_wt   => $dG_wt,
              dG_mt  => $dG_mt,
              ddG   => $ddG,
              dG_wt_sd  => $dG_wt_sd,
              dG_mt_sd   => $dG_mt_sd,
              ddG_sd  => $ddG_sd
            };

            $data->{evidence} = $evidence if (defined $evidence && $item eq "int");

            $formatted_output = $self->format_output($data, $item);
          }

          @$result_from_db{ keys %$formatted_output } = values %$formatted_output;
        }
      }
    }
  }

  return $result_from_db;
}

sub run {
  my ($self, $tva) = @_;
  
  my $result = {};

  # parse data file and generate result for tfbs
  if($self->{tfbs}){
    my $hash_from_file = $self->process_from_file($tva);
    @$result{ keys %$hash_from_file } = values %$hash_from_file;
  }
  # connect to db and generate result from matrix for the other type of analysis
  if($self->{motif} || $self->{int} || $self->{mod} || $self->{exp}){
    my $hash_from_db = $self->process_from_db($tva);
    @$result{ keys %$hash_from_db } = values %$hash_from_db;
  }

  return $self->{output_json} ? {"mutfunc" => $result} : $result;
}

sub parse_data {
  my ($self, $line) = @_;
  my ($c, $s, $ref, $alt, $del, undef, undef, undef, undef, undef, 
    $tf, $downstream, undef, undef, $g_strand, $s_strand, 
    $wt_score, $mt_score, undef, $ic_diff, undef, 
    $cells) = split /\t/, $line;
  
  return {
    chr => $c,
    start => $s,
    end => $s,
    ref => $ref,
    alt => $alt,
    result => {
      impact => $del, 
      tf => $tf, 
      downstream => $downstream, 
      g_strand => $g_strand,
      s_strand => $s_strand, 
      wt_score => $wt_score,
      mt_score => $mt_score,
      ic_diff => $ic_diff,
      cells => $cells
    }
  };
}

sub get_start { 
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;