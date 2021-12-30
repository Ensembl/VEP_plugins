=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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
 ./vep -i variations.vcf --plugin CADD,mutation_file=/full/path/to/intact/mutations.tsv,file=/full/path/to/genomic_coordinate/genomic_coordinates.vcf.gz
 ./vep -i variations.vcf --plugin CADD,mutation_file=/full/path/to/intact/mutations.tsv,file=/full/path/to/genomic_coordinate/genomic_coordinates.vcf.gz,minimal=1

=head1 DESCRIPTION

 A VEP plugin that retrieves molecular interaction data for variants from mutations
 data reported by IntAct database.
 
 Please cite the IntAct publication alongside the VEP if you use this resource:
 https://pubmed.ncbi.nlm.nih.gov/24234451/
 
 Pre-requisites:
 
 1) The IntAct mutation data file can be downloaded from -
 ftp://ftp.ebi.ac.uk/pub/databases/intact/current/various/mutations.tsv
 
 2) The genomic location for the HGVS ids needs to be generated to a file and needs 
 to be tabix-indexed. The current file format follows vcf format (see - 
 https://en.wikipedia.org/wiki/Variant_Call_Format) or needs to be tab-limited including 
 the following 5 fields in the order -
  CHR	POS	ID	REF	ALT
 
 One way to do this is using vep - 
  a) first extract the HGVS ids from the IntAct data file - 
  grep -v ^"#" mutations.tsv | cut -d'    ' -f2 > HGVS_in_mutations.txt

  b) second run vep on this HGVS ids (must be run with database) -
  ./vep --i HGVS_in_mutations.txt --o genomic_coordinates.vcf --vcf --database --db_version 104
  
  c) tabix index the generated vcf file -
  (grep ^"#" genomic_coordinates.vcf; grep -v ^"#" genomic_coordinates.vcf | sort -k1,1 -k2,2n) | bgzip > genomic_coordinates.vcf.gz
  tabix -p vcf genomic_coordinates.vcf.gz

 3) As you have already noticed, tabix utility must be installed in your path to use this plugin.
 
 Options are passed to the plugin as key=value pairs:

 file				: (mandatory) Path to tabix-indexed genomic coordinate file
 mutation_file			: (mandatory) Path to IntAct data file
 default			: (redundant) Set value to 1 to include feature_type, interaction_participants, 
				  pmid, and interaction_ac in the output.
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
 ap_organism			: Set value to 1 to include Affected protein organism in the output
 interaction_participants	: Set value to 1 to include Interaction participants in the output
 pmid				: Set value to 1 to include PubMedID in the output
 figure_legend			: Set value to 1 to include Figure legend in the output
 interaction_ac			: (redundant) Set value to 1 to include Interaction AC in the output. Adding
 				  this option redundant as it is already included by default or mimimal options.
 
 See what this options mean - https://www.ebi.ac.uk/intact/download/datasets#mutations

 Output:

 Output contains each of the selected field in order and separeted by comma. For example, a default output
 for variation O60701:p.Val132del can be -

 IntAct=mutation with no effect(MI:2226),uniprotkb:O60701(protein(MI:0326), 9606 - Homo sapiens),32001716,EBI-25431855

 Notice:
 
 The plugin is tested on ... 
 
=cut

package IntAct;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

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

  $self->add_file($param_hash->{file});
  
  $self->{interaction_ac} = 1;
  
  $self->{feature_type} = 1 unless defined $param_hash->{minimal};
  $self->{pmid} = 1 unless defined $param_hash->{minimal};
  $self->{interaction_participants} = 1 unless defined $param_hash->{minimal};
   
  $self->{feature_short_label} = 1 if defined $param_hash->{feature_short_label};
  $self->{feature_range} = 1 if defined $param_hash->{feature_range};
  $self->{original_sequence} = 1 if defined $param_hash->{original_sequence};
  $self->{resulting_sequence} = 1 if defined $param_hash->{resulting_sequence};
  $self->{feature_type} = 1 if defined $param_hash->{feature_type};
  $self->{feature_annotation} = 1 if defined $param_hash->{feature_annotation};
  $self->{ap_ac} = 1 if defined $param_hash->{ap_ac};
  $self->{ap_symbol} = 1 if defined $param_hash->{ap_symbol};
  $self->{ap_full_name} = 1 if defined $param_hash->{ap_full_name};
  $self->{ap_organism} = 1 if defined $param_hash->{ap_organism};
  $self->{figure_legend} = 1 if defined $param_hash->{figure_legend};

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;
  
  my %header;

  $header{"IntAct"} = "Molecular interaction data from IntAct database. Output is separated by , and includes: ";

  $header{"IntAct"} .= "Feature AC- Accession number for that particular mutation feature, " if $self->{feature_ac};
  $header{"IntAct"} .= "Feature short label- Human-readable short label summarizing the amino acid changes and their positions (HGVS compliant), " if $self->{feature_short_label};
  $header{"IntAct"} .= "Feature range(s)- Position(s) in the protein sequence affected by the mutation, " if $self->{feature_range};
  $header{"IntAct"} .= "Original sequence- Wild type amino acid residue(s) affected in one letter code, " if $self->{original_sequence};
  $header{"IntAct"} .= "Resulting sequence- Replacement sequence (or deletion) in one letter code, " if $self->{resulting_sequence};
  $header{"IntAct"} .= "Feature type- Mutation type following the PSI-MI controlled vocabularies, " if $self->{feature_type};
  $header{"IntAct"} .= "Feature annotation- Specific comments about the feature that can be of interest, " if $self->{feature_annotation};
  $header{"IntAct"} .= "Affected protein AC- Affected protein identifier (preferably UniProtKB accession, if available), " if $self->{ap_ac};
  $header{"IntAct"} .= "Affected protein symbol- As given by UniProtKB, " if $self->{ap_symbol};
  $header{"IntAct"} .= "Affected protein full name- As given by UniProtKB, " if $self->{ap_full_name};
  $header{"IntAct"} .= "Affected protein organism- TaxID and species name as given by UniProtKB, " if $self->{ap_organism};
  $header{"IntAct"} .= "Interaction participants- Identifiers for all participants in the affected interaction along with their species and molecule type between brackets, " if $self->{interaction_participants};
  $header{"IntAct"} .= "PubMedID- Reference to the publication where the interaction evidence was reported, " if $self->{pmid};
  $header{"IntAct"} .= "Figure legend- Reference to the specific figures in the paper where the interaction evidence was reported, " if $self->{figure_legend};
  $header{"IntAct"} .= "Interaction AC- Interaction accession within IntAct databases, " if $self->{interaction_ac};
  
  return \%header;
}

# match lines of IntAct data file using HGVS id
sub _match_id {
  my ($self, $id) = @_;

  my $param_hash = $self->params_to_hash();
  my $mutation_file = $param_hash->{mutation_file};

  my @matches;

  open my $f, $mutation_file or warning("WARNING: couldn't open file $mutation_file\n");
  while(<$f>) {
    if(/$id/) {
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
    or warning("WARNING: cannot parse intact data.\n");
    
  return \%parsed_data;
}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;

  # get allele
  my $allele = $tva->variation_feature_seq;
  
  return {} unless $allele =~ /^[ACGT-]+$/;

  my @data =  @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})};
  
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
    
    my $intact_matches = $self->_match_id($_->{id});
    
    return {} if (scalar @{ $intact_matches } != 1);

    my $parsed_data = $self->_parse_intact_data($intact_matches->[0]);
   
    my @result;

    my $i = 0;
    while (defined $valid_fields->{$i}){
      my $field = $valid_fields->{$i++};
      if(defined $self->{$field}) {
        push @result, $parsed_data->{$field} or warning("WARNING: failed to push result.\n");      
      }
    }    
    
    return {"IntAct" => [@result]} if (@result);
  }

  return {};
}

sub parse_data {
  my ($self, $line) = @_;
  my ($c, $s, $id, $ref, $alt) = split /\t/, $line ;
  
  # do VCF-like coord adjustment for mismatched subs
  my $e = ($s + length($ref)) - 1;
  if(length($alt) != length($ref)) {
    my $first_ref = substr($ref, 0, 1);
    my $first_alt = substr($alt, 0, 1);
    if ($first_ref eq $first_alt) {
      $s++;
      $ref = substr($ref, 1);
      $alt = substr($alt, 1);
      $ref ||= '-';
      $alt ||= '-';
    }
  }
  
  return {
    ref => $ref,
    alt => $alt,
    start => $s,
    end => $e,
    id => $id
  };
}

sub get_start {
  print $_, "\n";
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
