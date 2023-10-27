=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

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

 Paralogues

=head1 SYNOPSIS

 mv Paralogues.pm ~/.vep/Plugins

 # Download Ensembl paralogue annotation (if not available in current directory)
 # and fetch variants from Ensembl API (requires database access)
 ./vep -i variations.vcf --plugin Paralogues

 # Download Ensembl paralogue annotation (if not available in current directory)
 # and fetch variants from custom VCF file
 ./vep -i variations.vcf --plugin Paralogues,vcf=/path/to/file.vcf

 # Fetch Ensembl variants in paralogue proteins using only the Ensembl API
 # (requires database access)
 ./vep -i variations.vcf --database --plugin Paralogues,mode=remote

=head1 DESCRIPTION

 A VEP plugin that fetches variants overlapping the genomic coordinates of amino
 acids aligned between paralogue proteins. This is useful to predict
 pathogenicity of variants in paralogue positions.

 This plugin automatically downloads paralogue annotation from Ensembl databases
 if not available in the current directory (by default). Use argument `dir` to
 change the directory of the paralogue annotation:
   --plugin Paralogues,dir=/path/to/dir

 It is also possible to point to a custom tabix-indexed TSV file by using
 argument `paralogues`:
   --plugin Paralogues,paralogues=file.tsv.gz
   --plugin Paralogues,dir=/path/to/dir,paralogues=paralogues_file.tsv.gz
   --plugin Paralogues,paralogues=/path/to/dir/paralogues_file.tsv.gz

 The overlapping variants can be fetched from Ensembl API (by default) or a
 custom VCF file. To fetch variants from a VCF, use argument `vcf`:
   --plugin Paralogues,vcf=/path/to/file.vcf

 To avoid downloading data locally, the plugin also has a remote mode that
 optionally supports a custom VCF file:
   --plugin Paralogues,mode=remote
   --plugin Paralogues,mode=remote,vcf=/path/to/file.vcf

 The tabix utility must be installed in your path to use this plugin.

=cut

package Paralogues;

use strict;
use warnings;

use Bio::SimpleAlign;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

## GENERATE PARALOGUE ANNOTATION -----------------------------------------------

sub _prepare_filename {
  my $self = shift;
  my $config = $self->{config};
  my $reg    = $config->{reg};

  # Prepare file name based on species, database version and assembly
  my $pkg      = __PACKAGE__.'.pm';
  my $species  = $config->{species};
  my $version  = $config->{db_version} || $reg->software_version;
  my @name     = ($pkg, $species, $version);

  if( $species eq 'homo_sapiens' || $species eq 'human'){
    my $assembly = $config->{assembly} || $config->{human_assembly};
    die "specify assembly using --assembly [assembly]\n" unless defined $assembly;
    push(@name, $assembly) if defined $assembly;
  }
  return join("_", @name) . ".tsv.gz";
}

sub _get_homology_adaptor {
  my $self   = shift;
  my $config = $self->{config};
  my $reg    = $config->{reg};
  return $self->{ha} if $self->{ha};

  # reconnect to DB without species param
  if ($config->{host}) {
    $reg->load_registry_from_db(
      -host       => $config->{host},
      -user       => $config->{user},
      -pass       => $config->{password},
      -port       => $config->{port},
      -db_version => $config->{db_version},
      -no_cache   => $config->{no_slice_cache},
    );
  }

  return $reg->get_adaptor( "multi", "compara", "homology" );
}

sub _get_method_link_species_set_id {
  my $ha = shift;
  my $species = shift;
  my $type = shift;

  # Get ID corresponding to species and paralogues
  my @query = qq/
    SELECT method_link_species_set_id
    FROM method_link_species_set
    JOIN method_link USING (method_link_id)
    JOIN species_set ss USING (species_set_id)
    JOIN genome_db gdb ON ss.genome_db_id = gdb.genome_db_id
    WHERE gdb.name = "$species" AND type = "$type";
  /;
  my $sth = $ha->db->dbc->prepare(@query, { mysql_use_result => 1 });
  $sth->execute();
  my $id = @{$sth->fetchrow_arrayref}[0];
  $sth->finish();
  return $id;
}

sub _write_paralogue_annotation {
  my ($sth, $file) = @_;

  # Open lock
  my $lock = "$file\.lock";
  open LOCK, ">$lock" or
    die "ERROR: $file not found and cannot write to lock file $lock\n";
  print LOCK "1\n";
  close LOCK;

  open OUT, " | bgzip -c > $file" or die "ERROR: cannot write to file $file\n";

  # write header
  my @header = qw(homology_id chr start end strand stable_id version
                  perc_cov perc_id perc_pos cigar_line
                  paralogue_chr paralogue_start paralogue_end paralogue_strand
                  paralogue_stable_id paralogue_version paralogue_cigar_line);
  print OUT "#", join("\t", @header), "\n";

  while (my $line = $sth->fetchrow_arrayref()) {
    print OUT join("\t", @$line), "\n";
  }

  close OUT;
  unlink($lock);
  return $file;
}

sub _generate_paralogue_annotation {
  my $self = shift;
  my $file = shift;

  my $config  = $self->{config};
  my $species = $config->{species};
  my $reg     = $config->{reg};

  die("ERROR: Cannot generate paralogue annotation in offline mode\n") if $config->{offline};
  die("ERROR: Cannot generate paralogue annotation in REST mode\n") if $config->{rest};

  print "### Paralogues plugin: Querying Ensembl compara database (this may take a few minutes)\n" unless $config->{quiet};

  $self->{ha} ||= $self->_get_homology_adaptor;
  my $mlss_id = _get_method_link_species_set_id($self->{ha}, $species, 'ENSEMBL_PARALOGUES');

  # Create paralogue annotaton
  my @query = qq/
    SELECT
      hm.homology_id,
      df.name AS chr,
      sm.dnafrag_start   AS start,
      sm.dnafrag_end     AS end,
      sm.dnafrag_strand  AS strand,
      sm.stable_id,
      sm.version,
      hm.perc_cov,
      hm.perc_id,
      hm.perc_pos,
      hm.cigar_line,
      df2.name           AS paralogue_chr,
      sm2.dnafrag_start  AS paralogue_start,
      sm2.dnafrag_end    AS paralogue_end,
      sm2.dnafrag_strand AS paralogue_strand,
      sm2.stable_id      AS paralogue_id,
      sm2.version        AS paralogue_version,
      hm2.cigar_line     AS paralogue_cigar_line

    -- Reference proteins
    FROM homology_member hm
    JOIN homology h USING (homology_id)
    JOIN seq_member sm USING (seq_member_id)
    JOIN dnafrag df USING (dnafrag_id)

    -- Paralogue proteins
    JOIN homology_member hm2 ON hm.homology_id = hm2.homology_id
    JOIN seq_member sm2 ON hm2.seq_member_id = sm2.seq_member_id
    JOIN dnafrag df2 ON sm2.dnafrag_id = df2.dnafrag_id

    WHERE method_link_species_set_id = ${mlss_id} AND
          sm.source_name = 'ENSEMBLPEP' AND
          sm2.stable_id != sm.stable_id
    ORDER BY df.name, sm.dnafrag_start, sm.dnafrag_end;
  /;
  my $sth = $self->{ha}->db->dbc->prepare(@query, { mysql_use_result => 1});
  $sth->execute();

  print "### Paralogues plugin: Writing to file\n" unless $config->{quiet};
  _write_paralogue_annotation($sth, $file);
  $sth->finish();

  print "### Paralogues plugin: Creating tabix index\n" unless $config->{quiet};
  system "tabix -s2 -b3 -e4 $file" and die "ERROR: tabix index creation failed\n";

  print "### Paralogues plugin: file ready at $file\n" unless $config->{quiet};
  return 1;
}

sub _get_database_homologies {
  my $self       = shift;
  my $transcript = shift;

  my $config  = $self->{config};
  my $species = $config->{species};
  my $reg     = $config->{reg};

  $self->{ha} ||= $self->_get_homology_adaptor;
  $self->{ga} ||= $reg->get_adaptor($species, 'core', 'gene');
  my $gene = $self->{ga}->fetch_by_stable_id( $transcript->{_gene_stable_id} );
  my $homologies = $self->{ha}->fetch_all_by_Gene(
    $gene, -METHOD_LINK_TYPE => 'ENSEMBL_PARALOGUES', -TARGET_SPECIES => $species);
  return $homologies;
}

## RETRIEVE PARALOGUES FROM ANNOTATION -----------------------------------------

sub _compose_alignment_from_cigar {
  my $seq = shift;
  my $cigar = shift;

  die "Unsupported characters found in CIGAR line: $cigar\n"
    unless $cigar =~ /^[0-9MD]+$/;

  my $aln_str = $seq;
  my $index = 0;
  while ($cigar =~ /(\d*)([MD])/g) {
    my $num = $1 || 1;
    my $letter = $2;
    substr($aln_str, $index, 0) = '-' x $num if $letter =~ /D/;
    $index += $num;
  }
  return $aln_str;
}

sub _create_SimpleAlign {
  my ($self, $homology_id, $protein, $ref_cigar, $paralogue, $par_cigar) = @_;

  my $par       = $paralogue->translation;
  my $par_id    = $par->stable_id;
  my $par_seq   = $par->seq;

  my $ref_id  = $protein->stable_id;
  my $ref_seq = $protein->seq;

  my $aln = Bio::SimpleAlign->new();
  $aln->id($homology_id);
  $aln->add_seq(Bio::LocatableSeq->new(
    -SEQ        => _compose_alignment_from_cigar($ref_seq, $ref_cigar),
    -ALPHABET   => 'protein',
    -START      => 1,
    -END        => length($ref_seq),
    -ID         => $ref_id,
    -STRAND     => 0
  ));
  $aln->add_seq(Bio::LocatableSeq->new(
    -SEQ        => _compose_alignment_from_cigar($par_seq, $par_cigar),
    -ALPHABET   => 'protein',
    -START      => 1,
    -END        => length($par_seq),
    -ID         => $par_id,
    -STRAND     => 0
  ));

  # add paralogue information to retrieve from cache later on
  my $key = $par_id . '_info';
  $aln->{$key}->{_chr}    = $paralogue->seq_region_name;
  $aln->{$key}->{_start}  = $par->genomic_start;
  $aln->{$key}->{_strand} = $paralogue->strand;

  return $aln;
}

sub _get_paralogues {
  my ($self, $vf, $translation) = @_;
  my $var_chr   = $vf->seq_region_name || $vf->{chr};
  my $var_start = $vf->start - 2;
  my $var_end   = $vf->end;

  # get translation info
  my $translation_id  = $translation->stable_id;
  my $translation_seq = $translation->seq;

  # get paralogues for this variant region
  my $file = $self->{paralogues};
  my $data = `tabix $file $var_chr:$var_start-$var_end`;
  my @data = split /\n/, $data;

  my $paralogues = [];
  for (@data) {
    my (
      $homology_id, $chr, $start, $end, $strand, $protein_id, $version,
      $perc_cov, $perc_id, $perc_pos, $cigar,
      $para_chr, $para_start, $para_end, $para_strand, $para_id, $para_version, 
      $para_cigar
    ) = split /\t/, $_;
    next unless $translation_id eq $protein_id;

    my $paralogue = $self->_get_transcript_from_translation(
      $para_id, $para_chr, $para_start, $para_end);
    my $aln = $self->_create_SimpleAlign($homology_id, $translation, $cigar,
                                         $paralogue, $para_cigar);
    push @$paralogues, $aln;
  }
  return $paralogues;
}

sub _get_paralogue_coords {
  my $self = shift;
  my $tva  = shift;
  my $aln  = shift;

  my $translation_id    = $tva->transcript->translation->stable_id;
  my $translation_start = $tva->base_variation_feature_overlap->translation_start;

  # avoid cases where $translation_start is located in stop codon
  return unless $translation_start <= $tva->transcript->translation->length;

  # identify paralogue protein
  my @proteins = keys %{ $aln->{'_start_end_lists'} };
  my $paralogue;
  if ($translation_id eq $proteins[0]) {
    $paralogue = $proteins[1];
  } elsif ($translation_id eq $proteins[1]) {
    $paralogue = $proteins[0];
  } else {
    return;
  }

  # get genomic coordinates for aligned residue in paralogue
  my $col    = $aln->column_from_residue_number($translation_id, $translation_start);
  my $coords = $aln->get_seq_by_id($paralogue)->location_from_column($col);
  return unless defined $coords and $coords->location_type eq 'EXACT';

  my $tr_info          = $aln->{$paralogue . '_info'};
  my $tr_chr           = $tr_info->{_chr}    if defined $tr_info;
  my $tr_genomic_start = $tr_info->{_start}  if defined $tr_info;
  my $tr_strand        = $tr_info->{_strand} if defined $tr_info;
  my $para_tr = $self->_get_transcript_from_translation(
    $paralogue, $tr_chr, $tr_genomic_start, $tr_strand);

  my ($para_coords) = $para_tr->pep2genomic($coords->start, $coords->start);
  my $chr   = $para_tr->seq_region_name;
  my $start = $para_coords->start;
  my $end   = $para_coords->end;
  return ($chr, $start, $end);
}

## GET DATA FROM ANNOTATION ----------------------------------------------------

sub _get_transcript_from_translation {
  my $self       = shift;
  my $protein_id = shift;
  my $chr        = shift;
  my $start      = shift;
  my $strand     = shift;

  die "No protein identifier given\n" unless defined $protein_id;
  return $self->{_cache}->{$protein_id} if $self->{_cache}->{$protein_id};
  my $config  = Bio::EnsEMBL::VEP::Config->new( $self->{config} );

  # try to get transcript from cache
  if (defined $chr and defined $start and defined $strand) {
    my $as = $self->{as};
    if (!defined $as) {
      # cache AnnotationSource
      my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({config => $config});
      $as = $self->{as} = $asa->get_all_from_cache->[0];
    }

    my (@regions, $seen, $min_max, $min, $max);
    my $cache_region_size = $as->{cache_region_size};
    my $up_down_size = defined $as->{up_down_size} ? $as->{up_down_size} : $as->up_down_size;
    ($chr, $min, $max, $seen, @regions) = $as->get_regions_from_coords(
      $chr, $start, $start, $min_max, $cache_region_size, $up_down_size, $seen);

    foreach my $region (@regions) {
      my ($chr, $range) = @$region;

      my $file = $as->get_dump_file_name(
        $chr,
        ($range * $cache_region_size) + 1,
        ($range + 1) * $cache_region_size
      );

      my @features = @{
        $as->deserialized_obj_to_features( $as->deserialize_from_file($file) )
      } if -e $file;

      for my $transcript (@features) {
        my $translation = $transcript->translation;
        if ($translation && $translation->stable_id eq $protein_id) {
          $self->{_cache}->{$protein_id} = $transcript;
          return $transcript;
        }
      }
    }
  }

  # get transcript from database if not returned yet
  if ($self->{config}->{offline}) {
    die "Translation $protein_id not cached; avoid using --offline to allow to connect to database\n";
  }

  my $species = $config->species;
  my $reg     = $config->registry;
  $self->{ta} ||= $reg->get_adaptor($species, 'core', 'translation');

  my $transcript = $self->{ta}->fetch_by_stable_id($protein_id)->transcript;
  $self->{_cache}->{$protein_id} = $transcript;
  return $transcript;
}

## PLUGIN ----------------------------------------------------------------------

sub new {
  my $class = shift;  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();

  my $params = $self->params_to_hash();
  my $config = $self->{config};

  # Check for custom VCF
  my $vcf = $params->{vcf};
  if (defined $vcf) {
    $self->{vcf} = $vcf;
    $self->add_file($vcf);
  } elsif ($config->{offline} || $config->{rest}) {
    my $mode = $config->{rest} ? "REST" : "offline";
    die("ERROR: Cannot fetch Ensembl variants in $mode mode; please define vcf argument in the Paralogues plugin\n");
  }

  # Check if paralogue annotation should be downloaded
  $self->{remote} = defined $params->{mode} && $params->{mode} eq 'remote';
  return $self if $self->{remote};

  # Generate paralogue annotation
  my $annot = defined $params->{paralogues} ? $params->{paralogues} : $self->_prepare_filename;

  my $dir = $params->{dir};
  if (defined $dir) {
    die "ERROR: directory $dir not found" unless -e -d $dir;
    $annot = File::Spec->catfile( $dir, $annot );
  }

  $self->_generate_paralogue_annotation($annot) unless (-e $annot || -e "$annot.lock");
  $self->{paralogues} = $annot;
  return $self;
}

sub feature_types {
  return ['Feature'];
}

sub get_header_info {
  my $self = shift;

  my $header = {
    PARALOGUE_VARIANT_ID          => 'Paralogue variant identifier',
    PARALOGUE_VARIANT_SEQNAME     => 'Paralogue variant sequence region name',
    PARALOGUE_VARIANT_START       => 'Paralogue variant start',
    PARALOGUE_VARIANT_ALLELES     => 'Paralogue variant alleles',
  };

  if ( $self->{vcf} ) {
    $header->{'PARALOGUE_VARIANT_INFO_*'} = 'Paralogue variant fields based on VCF INFO field';
  } else {
    $header->{'PARALOGUE_VARIANT_STRAND'}      = 'Paralogue variant strand';
    $header->{'PARALOGUE_VARIANT_END'}         = 'Paralogue variant end';
    $header->{'PARALOGUE_VARIANT_CLNSIG'}      = 'Paralogue variant clinical significance';
    $header->{'PARALOGUE_VARIANT_SOURCE'}      = 'Paralogue variant source';
    $header->{'PARALOGUE_VARIANT_CONSEQUENCE'} = 'Paralogue variant consequence type(s)';
  }
  return $header;
}

sub _join_results {
  my $self = shift;
  my $all_results = shift;
  my $res = shift;

  if ($self->{config}->{output_format} eq 'json' || $self->{config}->{rest}) {
    # Group results for JSON and REST
    $all_results->{"Paralogues"} = [] unless defined $all_results->{"Paralogues"};
    push(@{ $all_results->{"Paralogues"} }, $res);
  } else {
    # Create array of results per key
    for (keys %$res) {
      $all_results->{$_} = [] if !$all_results->{$_};
      push(@{ $all_results->{$_} }, $res->{$_} || "NA");
    }
  }
  return $all_results;
}

sub _summarise_vf {
  my $vf = shift;
  return {
    PARALOGUE_VARIANT_ID          => $vf->name,
    PARALOGUE_VARIANT_SEQNAME     => $vf->seq_region_name,
    PARALOGUE_VARIANT_STRAND      => $vf->strand,
    PARALOGUE_VARIANT_START       => $vf->seq_region_start,
    PARALOGUE_VARIANT_END         => $vf->seq_region_end,
    PARALOGUE_VARIANT_ALLELES     => $vf->allele_string,
    PARALOGUE_VARIANT_CLNSIG      => join('/', @{$vf->get_all_clinical_significance_states}),
    PARALOGUE_VARIANT_SOURCE      => $vf->source_name,
    PARALOGUE_VARIANT_CONSEQUENCE => $vf->display_consequence,
  };
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  my $translation_start = $tva->base_variation_feature_overlap->translation_start;
  return {} unless defined $translation_start;

  my $homologies = [];
  if ($self->{remote}) {
    $homologies = $self->_get_database_homologies($tva->transcript);
  } else {
    my $translation = $tva->transcript->translation;
    $homologies = $self->{_cache_homologies}->{$translation->stable_id} ||=
      $self->_get_paralogues($vf, $translation);
  }
  return {} unless @$homologies;

  my $all_results = {};
  for my $aln (@$homologies) {
    $aln = $aln->get_SimpleAlign if $aln->isa('Bio::EnsEMBL::Compara::Homology');
    next unless $aln->isa('Bio::SimpleAlign');

    my ($chr, $start, $end) = $self->_get_paralogue_coords($tva, $aln);
    next unless defined $chr and defined $start and defined $end;

    if (defined $self->{vcf}) {
      # get variants from custom VCF file
      my @data = @{$self->get_data($chr, $start - 2, $end)};
      $all_results = $self->_join_results($all_results, $_) for @data;
    } else {
      # get Ensembl variants from mapped genomic coordinates
      my $config  = $self->{config};
      my $reg     = $config->{reg};
      my $species = $config->{species};

      $self->{sa}  ||= $reg->get_adaptor($species, 'core', 'slice');
      $self->{vfa} ||= $reg->get_adaptor($species, 'variation', 'variationfeature');

      my $slice = $self->{sa}->fetch_by_region('chromosome', $chr, $start, $end);
      foreach my $var ( @{ $self->{vfa}->fetch_all_by_Slice($slice) } ) {
        my $res = _summarise_vf($var);
        $all_results = $self->_join_results($all_results, $res);
      }
    }
  }
  return $all_results;
}

sub parse_data {
  my ($self, $line, $file) = @_;

  # parse custom VCF data
  my ($chrom, $start, $id, $ref, $alt, $qual, $filter, $info) = split /\t/, $line;

  # VCF-like adjustment of mismatched substitutions for comparison with VEP
  if(length($alt) != length($ref)) {
    my $first_ref = substr($ref, 0, 1);
    my $first_alt = substr($alt, 0, 1);
    if ($first_ref eq $first_alt) {
      $start++;
      $ref = substr($ref, 1);
      $alt = substr($alt, 1);
      $ref ||= '-';
      $alt ||= '-';
    }
  }

  my $data = {
    PARALOGUE_VARIANT_SEQNAME => $chrom,
    PARALOGUE_VARIANT_START   => $start,
    PARALOGUE_VARIANT_ID      => $id,
    PARALOGUE_VARIANT_ALLELES => "$ref/$alt"
  };

  # include data from INFO fields
  for my $field ( split /;/, $info ) {
    my ($key, $value) = split /=/, $field;
    $data->{'PARALOGUE_VARIANT_INFO_' . $key} = $value;
  };
  return $data;
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
