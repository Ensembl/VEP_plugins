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

 GO

=head1 SYNOPSIS

 mv GO.pm ~/.vep/Plugins

 # automatically fetch GFF files with GO terms and annotate input with GO terms
 # not compatible with --offline option
 ./vep -i variations.vcf --plugin GO
 
 # set directory used to write and read GFF files with GO terms
 ./vep -i variations.vcf --plugin GO,dir=${HOME}/go_terms

 # annotate input with GO terms from custom GFF file
 ./vep -i variations.vcf --plugin GO,file=${HOME}/custom_go_terms.gff.gz

 # annotate input based on gene identifiers instead of transcripts/translations
 ./vep -i variations.vcf --plugin GO,match=gene

 # use remote connection (available for compatibility purposes)
 ./vep -i variations.vcf --plugin GO,remote

=head1 DESCRIPTION

 A VEP plugin that retrieves Gene Ontology (GO) terms associated with
 transcripts (e.g. GRCh38) or their translations (e.g. GRCh37) using custom
 GFF annotation containing GO terms.

 The custom GFF files are automatically created if the input file do not exist
 by querying the Ensembl core database, according to database version, species
 and assembly used in VEP. Note that automatic retrieval fails if using the
 --offline option.

 The GFF files containing the GO terms are saved to and loaded from the working
 directory by default. To change this, provide a directory path as an argument:
 
   --plugin GO,dir=${HOME}/go_terms

 If your GFF file has a custom name, please provide the filename directly:

   --plugin GO,file=${HOME}/custom_go_terms.gff.gz
 
 The GO terms can also be fetched by gene match (either gene Ensembl ID or
 gene symbol) instead:

   --plugin GO,match=gene
   --plugin GO,match=gene_symbol
 
 To create/use a custom GFF file, these programs must be installed in your path:
   * The GNU zgrep and GNU sort commands to create the GFF file.
   * The tabix and bgzip utilities to create and read the GFF file: check
     https://github.com/samtools/htslib.git for installation instructions.

 Alternatively, for compatibility purposes, the plugin allows to use a remote
 connection to the Ensembl API by using "remote" as a parameter. This method
 retrieves GO terms one by one at both the transcript and translation level.
 This is not compatible with any other parameters:

   --plugin GO,remote

=cut

package GO;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  my $config = $self->{config};
  my $reg = $config->{reg};
  $reg = 'Bio::EnsEMBL::Registry';
  
  # Check if parameter "remote" is provided to revert to old GO.pm functionality
  $self->{use_remote} = grep($_ eq "remote", @{$self->{params}});
  
  # Check if the tabix command is available
  if ( !$self->{use_remote} and !`which tabix 2>&1` =~ /tabix$/ ) {
    die "ERROR: command tabix not found in your path\n" if $self->{config}{offline};
    
    # Use remote connection if online and if the tabix command is not available
    warn "WARNING: command tabix not found in your path so 'remote' was enabled\n";
    $self->{use_remote} = 1;
  }
  
  die "ERROR: cannot run 'remote' in offline mode\n" if ( $self->{use_remote} and $self->{config}{offline} );
  
  if ( !$self->{use_remote} ) {
    # Read GO terms from GFF file -- based on Phenotypes.pm

    my ($file, $dir, $match);
    my $param_hash = $self->params_to_hash();
    if (%$param_hash) {
      $file  = $param_hash->{file};
      $dir   = $param_hash->{dir};
      $match = $param_hash->{match};

      my @modes = ('transcript', 'translation', 'gene', 'gene_symbol');
      die "match argument $match is not valid -- available options are: " .
        join(", ", @modes) . "\n" if defined $match and !grep(/^$match$/, @modes);
    } elsif ( @{ $self->{params} } ) {
      $dir   = $self->{params}->[0];
    }
    $match ||= "transcript";
    $self->{match} = $match;

    if (defined $dir) {
      $dir =~ s/\/?$/\//; # ensure path ends with slash
      die "ERROR: directory $dir not found\n" unless -e -d $dir;
    }
    $dir  ||= "";
    $file ||= $self->_prepare_filename($reg);

    # Create GFF file with GO terms from database if file does not exist
    $file = $dir . $file;
    $self->_generate_gff($file) unless (-e $file || -e $file.'.lock');

    print "### GO plugin: Retrieving GO terms from $file\n" unless $config->{quiet};
    $self->add_file($file);
    $self->get_user_params();
  } else {
    # Revert to old GO.pm functionality -- based on Conservation.pm
    print "### GO plugin: Retrieving GO terms from Ensembl API\n" unless $config->{quiet};
        
    if(!defined($self->{config}->{sa})) {
      my $species = $config->{species};
      $reg->load_registry_from_db(
        -host       => $config->{host},
        -user       => $config->{user},
        -pass       => $config->{password},
        -port       => $config->{port},
        -db_version => $config->{db_version},
        -species    => $species =~ /^[a-z]+\_[a-z]+/i ? $species : undef,
        -verbose    => $config->{verbose},
        -no_cache   => $config->{no_slice_cache},
      );
    }
  }
  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  my $description = "GO terms associated with transcript or protein product";
  if ($self->{use_remote}) {
    $description .= " ('remote' mode)";
  } elsif ($self->{match} =~ /gene/) {
    $description .= sprintf(" ('%s' mode)", $self->{match});
  }
  return { 'GO' => $description};
}

sub run {
  my ($self, $tva) = @_;
  if ($self->{use_remote}) {
    # Remote connection to database
    return $self->_remote_run($tva);
  } else {
    # Match data from GFF file
    my $tr            = $tva->transcript;
    my $seqname       = $tr->{slice}->{seq_region_name};
    my $start         = $tr->{start};
    my $end           = $tr->{end};

    my $id;
    if ($self->{match} eq 'gene') {
      $id = $tr->{_gene_stable_id};
    } elsif ($self->{match} eq 'gene_symbol') {
      $id = $tr->{_gene_symbol};
    } else {
      $id = $tr->{stable_id};
    }
    return {} unless defined $id;

    my @data = @{$self->get_data($seqname, $start, $end)};
    foreach (@data) {
      return $_->{result} if $_->{id} eq $id;
    }
  }
  return {};
}

sub parse_data {
  my ($self, $line) = @_;
  my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attributes) = split /\t/, $line;

  # Parse ID and GO terms from attributes column
  my $id = undef;
  my $go = undef;
  foreach my $pair(split /;/, $attributes) {
    my ($key, $value) = split /\=/, $pair;
    next unless defined($key) and defined($value);
    if ($key eq "ID") {
      $id = $value;
    } elsif ($key eq "Ontology_term") {
      $go = $value;
    }
  }

  my $res;
  my @go_terms = split(",", $go);
  if ($self->{config}->{output_format} eq 'json' || $self->{config}->{rest}) {
    # Group results for JSON and REST
    my @go_terms_json = ();
    for (@go_terms) {
      my @items = split(":", $_);
      push @go_terms_json, {
        go_term => $items[0] . ":" . $items[1],
        description => $items[2]
      };
    }
    $res = \@go_terms_json;
  } else {
    $res = \@go_terms;
  }
  
  return {
    seqname => $seqname,
    start => $start,
    end => $end,
    id => $id,
    result => {
      GO => $res
    }
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

sub _prepare_filename {
  my ($self, $reg) = @_;
  my $config = $self->{config};
  my $match  = $self->{match};

  # Prepare file name based on species, database version and assembly
  my $pkg      = __PACKAGE__.'.pm';
  my $species  = $config->{species};
  my $version  = $config->{db_version} || $reg->software_version;
  my @basename = ($pkg, $species, $version);

  if( $species eq 'homo_sapiens' || $species eq 'human'){
    my $assembly = $config->{assembly} || $config->{human_assembly};
    die "specify assembly using --assembly [assembly]\n" unless defined $assembly;
    push @basename, $assembly if defined $assembly;
  }
  push @basename, $match if $match =~ 'gene';
  return join("_", @basename) . ".gff.gz";
}

sub _generate_gff {
  my ($self, $file) = @_;
  my $match = $self->{match};

  my $config = $self->{config};
  die("ERROR: Cannot create GFF file in offline mode\n") if $config->{offline};
  die("ERROR: Cannot create GFF file in REST mode\n") if $config->{rest};
  
  # test bgzip
  die "ERROR: bgzip does not seem to be in your path\n" unless `which bgzip 2>&1` =~ /bgzip$/;

  print "### GO plugin: Creating $file from database\n" unless($config->{quiet});
  
  print "### GO plugin: Querying Ensembl core database ('$match' match)\n"
    unless $config->{quiet};
  my $species = $config->{species};
  my $ta = $self->{config}->{reg}->get_adaptor($species, 'Core', 'Transcript');
  die ("ERROR: Ensembl core database not available\n") unless defined $ta;

  # Check whether GO terms are set at transcript or translation level
  my $id = _get_GO_terms_id( $ta );
  my $join_translation_table = _starts_with($id, "translation") ?
    "JOIN translation ON translation.transcript_id = transcript.transcript_id" : "";

  my $feature;
  my $join_gene_table;
  my $select_col;
  if ($match eq 'gene') {
    $feature = 'gene';
    $select_col = 'gene.stable_id';
    $join_gene_table = "JOIN gene ON gene.gene_id = transcript.gene_id";
  } elsif ($match eq 'gene_symbol') {
    $feature = 'gene';
    $select_col = 'gs.display_label';

    my $hgnc_id = _get_HGNC_id( $ta );
    $join_gene_table = qq/
      JOIN gene ON gene.gene_id = transcript.gene_id
      JOIN (
        SELECT stable_id, display_label FROM gene
        JOIN object_xref ox2 ON gene.gene_id = ox2.ensembl_id
        JOIN xref x2 ON ox2.xref_id = x2.xref_id AND x2.external_db_id = $hgnc_id
      ) AS gs ON gs.stable_id = gene.stable_id
    /;
  } else {
    $feature = 'transcript';
    $select_col = 'transcript.stable_id';
    $join_gene_table = '';
  }
  my $Feature = ucfirst $feature;

  # Query database for each GO term and its description per transcript
  my @query = qq/
    SELECT DISTINCT
      sr.name AS seqname,
      REPLACE(db.db_name, " ", "_") AS source,
      "$Feature" AS feature,
      $feature.seq_region_start AS start,
      $feature.seq_region_end AS end,
      '.' AS score,
      IF($feature.seq_region_strand = 1, '+', '-') AS strand, 
      '.' AS frame,
      $select_col,
      x.display_label AS go_term,
      x.description AS go_term_description
      
    FROM transcript
    $join_translation_table
    JOIN object_xref ox ON $id = ox.ensembl_id
    JOIN xref x ON ox.xref_id = x.xref_id
    JOIN seq_region sr ON transcript.seq_region_id = sr.seq_region_id
    JOIN external_db db ON x.external_db_id = db.external_db_id
    $join_gene_table
    WHERE db.db_name = "GO" AND x.dbprimary_acc LIKE "GO:%"
    # the following order is a major downstream assumption
    ORDER BY $select_col, x.display_label
  /;
  my $sth = $ta->db->dbc->prepare(@query, { mysql_use_result => 1});
  $sth->execute();

  # Append all GO terms from the same feature and write to file
  print "### GO plugin: Writing to file\n" unless $config->{quiet};
  my $file_tmp = _write_GO_terms_to_file($sth, $file);
  $sth->finish();
  
  print "### GO plugin: Sorting file\n" unless $config->{quiet};
  system("(zgrep '^#' $file_tmp; LC_ALL=C zgrep -v '^#' $file_tmp | sort -k1,1 -k4,4n ) | bgzip -c > $file") and die("ERROR: sort failed\n");
  unlink($file_tmp);

  print "### GO plugin: Creating tabix index\n" unless $config->{quiet};
  system "tabix -p gff $file" and die "ERROR: tabix index creation failed\n";

  print "### GO plugin: GFF file ready!\n" unless $config->{quiet};
  return 1;
}

sub _starts_with {
  my ($string, $prefix) = @_;
  return rindex($string, $prefix, 0) == 0;
}

sub _get_HGNC_id {
  my ($ta) = @_;
  my @query = "SELECT external_db_id FROM external_db WHERE db_name = 'HGNC'";
  my $sth = $ta->db->dbc->prepare(@query, { mysql_use_result => 1});
  $sth->execute();
  my $id = @{$sth->fetchrow_arrayref}[0];
  $sth->finish();
  return $id;
}

sub _get_GO_terms_id {
  my ($ta) = @_;
  
  my @query = qq{
    SELECT ox.ensembl_object_type
    FROM ontology_xref go
    LEFT JOIN object_xref ox ON go.object_xref_id = ox.object_xref_id
    LIMIT 1;
  };
  my $sth = $ta->db->dbc->prepare(@query, { mysql_use_result => 1});
  $sth->execute();
  my $type = lc( @{$sth->fetchrow_arrayref}[0] );
  $sth->finish();
  return "$type.$type\_id";
}

sub _write_GO_terms_to_file {
  my ($sth, $file) = @_;
  my $file_tmp = $file.".tmp";
  
  # Open lock
  my $lock = "$file\.lock";
  open LOCK, ">$lock" or
    die "ERROR: $file not found and cannot write to lock file $lock\n";
  print LOCK "1\n";
  close LOCK;

  open OUT, " | bgzip -c > $file_tmp" or die "ERROR: cannot write to file $file_tmp\n";
  print OUT "##gff-version 1.10\n"; # GFF file header

  # For a single feature, append all of its GO terms to $info;
  # when there is a new feature, write $info to file and repeat
  my $info;
  my $previous_feature = "";
  while(my $line = $sth->fetchrow_arrayref()) {
    my $row = [ @$line ];
    my ($id, $go_term, $description) = splice(@$row, -3);

    if ($id ne $previous_feature) {
      # If not the same feature, write previous feature info to file
      print OUT $info."\n" if defined($info);

      # Set this new feture info
      $previous_feature = $id;
      $row = join("\t", map {defined($_) ? $_ : '.'} @$row);
      $info = $row."\tID=$id;Ontology_term=";
    } else {
      # Append comma before appending another GO term
      $info .= ","
    }

    if ( defined($description) ) {
      $description =~ s/ /_/g; # Replace spaces with underscores

      $description =~ s/,_/_-_/g; # Avoid commas followed by an underscore, e.g.:
      # GO:0045892:negative_regulation_of_transcription,_DNA-templated

      $description =~ s/,/-/g; # Avoid commas in other situtations, e.g.:
      # GO:0016316:phosphatidylinositol-3,4-bisphosphate_4-phosphatase_activity
    } else {
      $description = "";
    }

    # Append GO term and its description
    $info .= "$go_term:$description";
  }
  # Write info of last feature to file
  print OUT $info."\n" if defined($info);

  close OUT;
  unlink($lock);
  return $file_tmp;
}

sub _remote_run {
  my ($self, $tva) = @_;
  
  my $tr = $tva->transcript;
  return {} unless defined($tr);
  
  # Get GO terms at transcript and translation levels
  my $entries = $tr->get_all_DBLinks('GO');

  # Format ID and description of GO terms (and ignore duplicates)
  my @go_terms = _uniq( map {$_->display_id.':'.$_->description} @$entries );
  my $string = join(",", @go_terms);
  $string =~ s/\s+/\_/g;
  
  # Avoid returning empty GO terms
  return $string eq "" ? {} : { GO => $string };
}

sub _uniq {
  my %seen;
  grep !$seen{$_}++, @_;
}

1;
