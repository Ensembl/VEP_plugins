
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

 ProtVar

=head1 SYNOPSIS

 mv ProtVar.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin ProtVar,db=/FULL_PATH_TO/ProtVar.db
 ./vep -i variations.vcf --plugin ProtVar,db=/FULL_PATH_TO/ProtVar_data.db,stability=1

=head1 DESCRIPTION

 An Ensembl VEP plugin that retrieves data from ProtVar resource providing contexualised information for
 missense variation such destabilization of protein structure, overlapping protein pocket, 
 and protein-protein interaction interface.

 Please cite the ProtVar publication alongside Ensembl VEP if you use this resource:
 https://academic.oup.com/nar/article/52/W1/W140/7676839

 Pre-requisites:

 1) The data file. ProtVar SQLite db can be downloaded from - 
 https://ftp.ensembl.org/pub/current_variation/ProtVar/ProtVar_data.db

 2) If you are using --offline please provide a FASTA file as this plugin requires the
 translation sequence to function.

 Options are passed to the plugin as key=value pairs:

 db         : (mandatory) Path to SQLite database containing the data.
 stability  : Select this option to have protein stability related output. Output contains -
                ddG - Free energy change upon mutation
 pocket     : Select this option to have overlapping protein pocket related output. Output contains -
                id - Pocket id;
                score - Combined score measuring confidence in pocket (score > 800 high confidence and score > 900 very high confidence)
                MpLDDT - Mean pLDDT score of all the residues from AlphaFold2 model used to form the pocket
                energy - Pocket energy per volume (in kcal/mol unit)
                burriedness - Pocket burriedness (0.0 means entirely exposed and 1.0 means entierly burried)
                RoG - Pocket compactness in terms of radius of gyration  (in Angstrom unit)
                residues - Position of residues that form the pocket according to UniProt canonical (and AlphaFold) numbering.
 int      : Select this option to have protein-protein interface related output. Output contains -
                protein - Uniprot id of the protein with which the overlapping protein creates an interface
                pDockQ - A docking score for the interface, higher scores imply a more stable interface.

 By default all of the three type of outputs (stability, pocket, int) are provided. But if you want to have some selected type and not 
 all of them just select the relevant options.

=cut

package ProtVar;

use strict;
use warnings;
use DBI;
use Compress::Zlib;
use Digest::MD5     qw(md5_hex);
use List::MoreUtils qw(first_index);

use Bio::EnsEMBL::Utils::Exception           qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my @ALL_AAS = qw(A C D E F G H I K L M N P Q R S T V W Y);

my $field_order = {
    stability => [ "ddG" ],
    pocket => [ "id", "score", "MpLDDT", "energy", "burriedness", "RoG", "residues" ],
    int => [ "protein", "pDockQ" ]
};

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);

    my $param_hash = $self->params_to_hash();

    if ( $self->{config}->{offline} && !$self->{config}->{fasta} ) {
        die "The ProtVar plugin cannot annotate without translation sequence, please provide a FASTA if --offline is used\n";
    }

    # default behavior is to output all field
    $param_hash->{all} = 1
      if ( ( !defined $param_hash->{stability} )
        && ( !defined $param_hash->{pocket} )
        && ( !defined $param_hash->{int} ) );

    $self->{stability} = 1 if $param_hash->{stability} || $param_hash->{all};
    $self->{pocket}    = 1 if $param_hash->{pocket}    || $param_hash->{all};
    $self->{int}       = 1 if $param_hash->{int}       || $param_hash->{all};

    die "ERROR: please provide the SQLite database using 'db' parameter\n" if ( !( defined $param_hash->{db} ) );
    $self->{db} = $param_hash->{db};

    if ($self->{config}->{output_format} && $self->{config}->{output_format} eq "json") {
        $self->{output_json} = 1;
    }
    if ( $self->{config}->{rest} ) {
        $self->{output_json} = 1;
    }

    $self->{initial_pid} = $$;

    $self->{dbh} = DBI->connect( "dbi:SQLite:dbname=" . $self->{db}, "", "" );
    $self->{get_sth} = $self->{dbh}->prepare("SELECT md5, item, matrix FROM consequences WHERE md5 = ?");

    return $self;
}

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    my ($self) = shift;

    my %header;

    if ( defined $self->{stability} ) {
        $header{ProtVar_stability} = "Impact on protein stability in ddG enegry difference.";
    }

    if ( defined $self->{pocket} ) {
        $header{ProtVar_pocket} = "Information about overlapping protein pocket. Output field(s) include: ";
        $header{ProtVar_pocket} .=
          $self->{config}->{output_format} eq "vcf"
          ? "(fields are separated by '&') "
          : "(fields are separated by ',') ";
        $header{ProtVar_pocket} .= "id - Pocket id, ";
        $header{ProtVar_pocket} .= "score - Combined score measuring confidence in pocket (score > 800 high confidence and score > 900 very high confidence), ";
        $header{ProtVar_pocket} .= "MpLDDT - Mean pLDDT score of all the residues from AlphaFold2 model used to form the pocket, ";
        $header{ProtVar_pocket} .= "energy - Pocket energy per volume (in kcal/mol unit), ";
        $header{ProtVar_pocket} .= "burriedness - Pocket burriedness (0.0 means entirely exposed and 1.0 means entierly burried), ";
        $header{ProtVar_pocket} .= "RoG - Pocket compactness in terms of radius of gyration  (in Angstrom unit), ";
        $header{ProtVar_pocket} .= "residues - Position of residues that form the pocket according to UniProt canonical (and AlphaFold) numbering.";
    }

    if ( defined $self->{int} ) {
        $header{ProtVar_int} = "Information about overlapping protein-protein interface . Output field(s) include: ";
        $header{ProtVar_int} .=
          $self->{config}->{output_format} eq "vcf"
          ? "(fields are separated by '&') "
          : "(fields are separated by ',') ";
        $header{ProtVar_int} .= "protein_id - ID of the protein that overlapping protien interface with, ";
        $header{ProtVar_int} .= "pDockQ - A docking score for the interface, higher scores imply a more stable interface.";
    }

    return \%header;
}

sub expand_matrix {
    my ($matrix) = @_;
    my $expanded_matrix = Compress::Zlib::memGunzip($matrix)
      or throw("Failed to gunzip: $gzerrno");

    return $expanded_matrix;
}

sub retrieve_item_value {
    my ( $matrix, $pos, $aa, $tot_packed_len, $item ) = @_;

    my $item_value;
    if ( $item eq 'stability' ) {
        $item_value = substr $matrix,
          $pos * $tot_packed_len * 20 + $aa * $tot_packed_len, $tot_packed_len;
    }
    if ( $item eq 'pocket' || $item eq 'int' ) {
        $item_value = substr $matrix, $pos * $tot_packed_len, $tot_packed_len;
    }

    return $item_value;
}

sub parse_stability {
    my ($item_value) = @_;

    my $ddG = unpack "A8", $item_value;

    $ddG = undef if $ddG eq "10000000";

    return $ddG;
}

sub parse_pocket {
    my ($item_value) = @_;

    my $id = unpack "v", $item_value;

    $id = undef if $id == 0xFFFF;

    return $id;
}

sub parse_int {
    my ($item_value) = @_;

    my ( $protein, $pDockQ ) = unpack "A12A8", $item_value;

    $protein = undef if $protein eq "undefined";
    $pDockQ  = undef if $pDockQ eq "10000000";

    return $protein, $pDockQ;
}

sub format_output {
    my ( $self, $data, $item ) = @_;

    my $result = {};

    if ( $self->{output_json} ) {
        my %hash;
        %hash = map { $_ => $data->{$_} } keys %$data;
        $result->{$item} = \%hash;
    }
    else {
        my $key = "ProtVar_" . $item;
        
        # replace comma in pocket residue position to avoid delimiter clash
        if ($data->{residues}) {
            $data->{residues} =~ s/,/p/g;
            $data->{residues} = 'p'.$data->{residues};
        }
        
        $result->{$key} = join(",", map { $data->{$_} } @{ $field_order->{$item} });
    }

    return $result;
}

sub process_from_db {
    my ( $self, $tva ) = @_;

    # get the trascript related to the variant
    my $tr = $tva->transcript;

    # get the translation
    my $translation = $tr->translate;
    return {} unless $translation;

    # get the md5 hash of the peptide sequence
    my $md5 = md5_hex( $translation->seq );

    # forked, reconnect to DB
    if ( $$ != $self->{initial_pid} ) {
        $self->{dbh} = DBI->connect( "dbi:SQLite:dbname=" . $self->{db}, "", "" );
        $self->{get_sth} = $self->{dbh}->prepare("SELECT md5, item, matrix FROM consequences WHERE md5 = ?");

        # set this so only do once per fork
        $self->{initial_pid} = $$;
    }
    $self->{get_sth}->execute($md5);

    my $result_from_db = {};
    while ( my $arrayref = $self->{get_sth}->fetchrow_arrayref ) {
        my $item   = $arrayref->[1];
        my $matrix = $arrayref->[2];

        # expand the compressed matrix
        my $expanded_matrix = expand_matrix($matrix);

        # position and peptide to retrieve value from matrix
        my $pos     = $tva->transcript_variation->translation_start;
        my $peptide = $tva->peptide;

        next if ( $item eq "stability" && !defined $peptide );

        # in matrix position is 0 indexed
        $pos--;

        # we need position of peptide in the ALL_AAS array
        my $peptide_number = ( first_index { $_ eq $peptide } @ALL_AAS );

        # get matrix for each item value and parse them
        # stability
        if ( $item eq "stability" && $self->{stability} ) {

            # get the item value for specific pos and amino acid from matrix
            my $tot_packed_len = 8;
            my $item_value = retrieve_item_value( $expanded_matrix, $pos, $peptide_number, $tot_packed_len, 'stability' );

            if ($item_value) {

                # parse the item value from matrix
                my $ddG = parse_stability($item_value);

                # format the output
                if ( $ddG ) {
                    my $data = { 
                        ddG => $ddG 
                    };

                    my $formatted_output = $self->format_output( $data, $item );
                    @$result_from_db{ keys %$formatted_output } = values %$formatted_output;
                }
            }
        }

        # pocket
        elsif ( $item eq "pocket" && $self->{pocket} ) {

            # get the item value for specific pos and amino acid from matrix
            my $tot_packed_len = 2;
            my $item_value = retrieve_item_value( $expanded_matrix, $pos, $peptide_number, $tot_packed_len, 'pocket' );

            if ($item_value) {

                # parse the item value from matrix
                my $id = parse_pocket($item_value);

                # format the output
                if ($id) {
                    my $sth2 = $self->{dbh}->prepare( "SELECT score, MpLDDT, energy, burriedness, RoG, residues FROM pocket_attribs WHERE md5 = ? AND pocket_id = ?" );
                    $sth2->execute($md5, $id);
                    my ( $score, $MpLDDT, $energy, $burriedness, $RoG, $residues ) = @{ $sth2->fetchrow_arrayref };
                    
                    my $data = {
                        id => 'P' . $id,
                        score => $score,
                        MpLDDT => $MpLDDT,
                        energy => $energy,
                        burriedness => $burriedness,
                        RoG => $RoG,
                        residues => $residues
                    };

                    my $formatted_output = $self->format_output( $data, $item );
                    @$result_from_db{ keys %$formatted_output } = values %$formatted_output;
                }
            }
        }

        # int
        elsif ( $item eq "int" && $self->{int} ) {

            # get the item value for specific pos and amino acid from matrix
            my $tot_packed_len = 20;
            my $item_value = retrieve_item_value( $expanded_matrix, $pos, $peptide_number, $tot_packed_len, 'int' );

            if ($item_value) {

                # parse the item value from matrix
                my ( $protein, $pDockQ ) = parse_int($item_value);

                # format the output
                if ( $protein && $pDockQ ) {
                    my $data = {
                        protein => $protein,
                        pDockQ  => $pDockQ
                    };

                    my $formatted_output = $self->format_output( $data, $item );
                    @$result_from_db{ keys %$formatted_output } = values %$formatted_output;
                }
            }
        }
    }

    return $result_from_db;
}

sub run {
    my ( $self, $tva ) = @_;

    my $result = {};

    my $hash_from_db = $self->process_from_db($tva);
    @$result{ keys %$hash_from_db } = values %$hash_from_db;

    return {} unless %$result;
    return $self->{output_json} ? { "ProtVar" => $result } : $result;
}

1;
