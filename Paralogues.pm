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

 Paralogues

=head1 SYNOPSIS

 mv Paralogues.pm ~/.vep/Plugins

 # Find paralogue regions of all input variants using Ensembl paralogue annotation
 # (automatically created if not in current directory) and fetch variants within
 # those regions from VEP cache and whose clinical significance partially
 # matches 'pathogenic'
 ./vep -i variations.vcf --cache --plugin Paralogues

 # Find paralogue regions of input variants using Ensembl paralogue annotation
 # (automatically created if not in current directory) and fetch variants within
 # those regions from a custom VCF file (regardless of their clinical significance)
 ./vep -i variations.vcf --cache --plugin Paralogues,vcf=/path/to/file.vcf,clnsig=ignore

 # Same using a custom VCF file but filtering for 'pathogenic' variants
 ./vep -i variations.vcf --cache --plugin Paralogues,vcf=/path/to/file.vcf,clnsig_col=CLNSIG

 # Fetch all Ensembl variants in paralogue proteins using only the Ensembl API
 # (requires database access)
 ./vep -i variations.vcf --database --plugin Paralogues,mode=remote,clnsig='ignore'

=head1 DESCRIPTION

 A VEP plugin that fetches variants overlapping the genomic coordinates of amino
 acids aligned between paralogue proteins. This is useful to predict the
 pathogenicity of variants in paralogue positions.

 This plugin determines paralogue regions for a variant based on Ensembl
 paralogue annotation. These versatile annotations allow to look up paralogue
 regions for all variants from any species with Ensembl paralogues.

 After retrieving the paralogue regions, this plugin fetches variants
 overlapping those regions from one of the following sources (by this order):
   1. Custom VCF via the 'vcf' parameter
   2. VEP cache (in cache/offline mode)
   3. Ensembl API (in database mode)

 This plugin automatically downloads paralogue annotation from Ensembl databases
 if not available in the current directory (by default). Use option `dir` to
 change the directory of the paralogue annotation:
   --plugin Paralogues
   --plugin Paralogues,dir=/path/to/dir

 It is also possible to point to a custom tabix-indexed TSV file by using
 option `paralogues`:
   --plugin Paralogues,paralogues=file.tsv.gz
   --plugin Paralogues,dir=/path/to/dir,paralogues=paralogues_file.tsv.gz
   --plugin Paralogues,paralogues=/path/to/dir/paralogues_file.tsv.gz

 Returned variants can be filtered based on clinical significance by using
 option `clnsig` (use keyword 'ignore' to avoid this filtering):
   --plugin Paralogues,clnsig=ignore
   --plugin Paralogues,clnsig=pathogenic,clnsig_match=partial
   --plugin Paralogues,clnsig='likely pathogenic',clnsig_match=exact
   --plugin Paralogues,vcf=/path/to/file.vcf.gz,clnsig=benign,clnsig_col=CLNSIG

 The overlapping variants can be fetched from Ensembl API (by default) or a
 custom tabix-indexed VCF file. You can point to a VCF file using option `vcf`
 and input a colon-delimited list of INFO fields in `fields`:
   --plugin Paralogues,vcf=/path/to/file.vcf.gz,clnsig_col=CLNSIG
   --plugin Paralogues,vcf=/path/to/file.vcf.gz,clnsig=ignore,fields=identifier:alleles:CLNSIG:CLNVI:GENEINFO

 Options are passed to the plugin as key=value pairs:
   dir          : Directory with paralogue annotation (the annotation is created
                  in this folder if the required files do not exist)
   paralogues   : Tabix-indexed TSV file with paralogue annotation (the
                  annotation is created with this filename, if the file does not
                  exist); if set to 'remote', paralogue annotation is fetched
                  from Ensembl API
   min_perc_cov : Minimum alignment percentage of the peptide associated with
                  the input variant (default: 0)
   min_perc_pos : Minimum percentage of positivity (similarity) between both
                  homologues (default: 50)
   regions      : Boolean value to return regions used to look up paralogue
                  variants (default: 1)

   vcf          : Custom Tabix-indexed VCF file with information for overlapping 
                  variants (by default, overlapping variants are fetched from
                  VEP cache in cache/offline mode or Ensembl API in database mode)
   fields       : Colon-separated list of information from paralogue variants to
                  output (default: 'identifier:alleles:clinical_significance');
                  keyword 'all' can be used to print all fields; available
                  fields include 'identifier', 'chromosome', 'start', 'alleles',
                  'perc_cov', 'perc_pos', and 'clinical_significance' (if
                  `clnsig_col` is defined for custom VCF); additional fields
                  are available depending on variant source:
                    - VEP cache: 'end' and 'strand'
                    - Ensembl API: 'end', 'strand', 'source', 'consequence' and
                      'gene_symbol'
                    - Custom VCF: 'quality', 'filter' and name of INFO fields

   clnsig       : Clinical significance term to filter variants (default:
                  'pathogenic'); use 'ignore' to fetch all paralogue variants,
                  regardless of clinical significance
   clnsig_match : Type of match when filtering variants based on option
                  `clnsig`: 'partial' (default), 'exact' or 'regex'
   clnsig_col   : Column name containing clinical significance in custom VCF
                  (required with `vcf` option and if `clnsig` is not 'ignore')

 The tabix utility must be installed in your path to read the paralogue
 annotation and the custom VCF file.

=cut

package Paralogues;

use strict;
use warnings;

use List::Util qw(any);
use File::Basename;
use Compress::Zlib;

use Bio::SimpleAlign;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::IO::Parser::VCF4Tabix;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

our @VCF_FIELDS = (
  'identifier',
  'chromosome',
  'start',
  'alleles',
  'quality',
  'filter',
  'clinical_significance',
  'perc_cov',
  'perc_pos',
);

our @CACHE_FIELDS = (
  'identifier',
  'chromosome',
  'start',
  'end',
  'strand',
  'alleles',
  'clinical_significance',
  'perc_cov',
  'perc_pos',
);

our @API_FIELDS = (
  @CACHE_FIELDS,
  'source',
  'consequence',
  'gene_symbol',
);

## GENERATE PARALOGUE ANNOTATION -----------------------------------------------

sub _prepare_filename {
  my $self = shift;
  my $config = $self->{config};

  # Prepare file name based on species, database version and assembly
  my $pkg      = __PACKAGE__.'.pm';
  my $species  = $config->{species};
  my $version  = $config->{db_version} || 'Bio::EnsEMBL::Registry'->software_version;
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

  if (!defined $self->{ha}) {
    $self->{ha} = $reg->get_adaptor( "multi", "compara", "homology" );
  }

  if (!defined $self->{ha}) {
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
    $self->{ha} = $reg->get_adaptor( "multi", "compara", "homology" );
  }
  return $self->{ha};
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

  my $mlss_id = _get_method_link_species_set_id($self->_get_homology_adaptor, $species, 'ENSEMBL_PARALOGUES');

  # Create paralogue annotation
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

  $self->{ga} ||= $reg->get_adaptor($species, 'core', 'gene');
  my $gene = $self->{ga}->fetch_by_stable_id( $transcript->{_gene_stable_id} );
  my $homologies = $self->_get_homology_adaptor->fetch_all_by_Gene(
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
  my ($self, $homology_id, $protein, $ref_cigar, $paralogue, $par_cigar,
      $perc_cov, $perc_pos) = @_;

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

  # add alignment stats
  $aln->{_stats}->{ref_perc_cov} = $perc_cov;
  $aln->{_stats}->{ref_perc_pos} = $perc_pos;

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
  my @data = @{$self->get_data($var_chr, $var_start, $var_end, $file)};

  my $paralogues = [];
  for (@data) {
    my (
      $homology_id, $chr, $start, $end, $strand, $protein_id, $version,
      $perc_cov, $perc_id, $perc_pos, $cigar,
      $para_chr, $para_start, $para_end, $para_strand, $para_id, $para_version, 
      $para_cigar
    ) = split /\t/, $_;
    next unless $translation_id eq $protein_id;
    next unless $perc_cov >= $self->{min_perc_cov};
    next unless $perc_pos >= $self->{min_perc_pos};

    my $paralogue = $self->_get_transcript_from_translation(
      $para_id, $para_chr, $para_start, $para_end);
    my $aln = $self->_create_SimpleAlign($homology_id, $translation, $cigar,
                                         $paralogue, $para_cigar,
                                         $perc_cov, $perc_pos);
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

sub _get_config {
  my $self = shift;
  if (!defined $self->{config_obj}) {
    $self->{config_obj} = Bio::EnsEMBL::VEP::Config->new( $self->{config} );
  }
  return $self->{config_obj};
}

sub _get_AnnotationSource {
  my $self = shift;
  my $filter = shift;

  my %match = (
    'Transcript'     => 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript',
    'Variation'      => 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation',
    'VariationTabix' => 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix',
  );

  if (!defined $self->{asa}) {
    # Cache all annotation sources
    my $cfg = $self->_get_config;
    $cfg->{_params}->{check_existing} = 1; # enable to fetch variants from VEP cache
    $self->{asa} = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({config => $cfg});
  }

  my $asa = $self->{asa};
  if (defined $asa && $asa->can('get_all_from_cache')) {
    for (@{$self->{asa}->get_all_from_cache}) {
      return $_ if ref $_ eq $match{$filter};
    }
  }
  return undef;
}

sub _get_transcript_from_translation {
  my $self       = shift;
  my $protein_id = shift;
  my $chr        = shift;
  my $start      = shift;
  my $strand     = shift;

  die "No protein identifier given\n" unless defined $protein_id;
  return $self->{_cache}->{$protein_id} if $self->{_cache}->{$protein_id};

  # try to get transcript from cache if enabled
  if ($self->{config}->{cache} && defined $chr && defined $start && defined $strand) {
    my $as = $self->_get_AnnotationSource('Transcript');
    die "ERROR: could not get transcripts from VEP cache" unless defined $as;

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

  my $config  = $self->_get_config;
  my $species = $config->species;
  my $reg     = $config->registry;
  $self->{ta} ||= $reg->get_adaptor($species, 'core', 'translation');

  my $transcript = $self->{ta}->fetch_by_stable_id($protein_id)->transcript;
  $self->{_cache}->{$protein_id} = $transcript;
  return $transcript;
}

## GET PARALOGUE VARIANTS BASED ON MATCHES OR ANNOTATION -----------------------

sub _get_paralogue_vars_from_annotation {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;

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
    my ($perc_cov, $perc_pos);
    if ($aln->isa('Bio::EnsEMBL::Compara::Homology')) {
      my $ref = $aln->get_all_Members->[0];
      $perc_cov = $ref->perc_cov;
      $perc_pos = $ref->perc_pos;
      $aln = $aln->get_SimpleAlign;
    } elsif ($aln->isa('Bio::SimpleAlign')) {
      $perc_cov = $aln->{_stats}->{ref_perc_cov};
      $perc_pos = $aln->{_stats}->{ref_perc_pos};
    } else {
      next;
    }

    my ($chr, $start, $end) = $self->_get_paralogue_coords($tva, $aln);
    next unless defined $chr and defined $start and defined $end;
    $all_results = $self->_get_paralogue_vars($chr, $start, $end, $perc_cov, $perc_pos, $all_results);
  }
  return $all_results;
}

## FETCH VARIANTS FROM SOURCE --------------------------------------------------

sub _fetch_cache_vars {
  # fetch variants from VEP cache
  my ($self, $chr, $start, $end) = @_;

  my ($as, $variants);
  if (defined($as = $self->_get_AnnotationSource('VariationTabix'))) {
    # code based on AnnotationSource::Cache::VariationTabix
    my $source_chr = $as->get_source_chr_name($chr);
    my $tabix_obj = $as->_get_tabix_obj($source_chr);
    next unless $tabix_obj;

    my $iter = $tabix_obj->query(sprintf("%s:%i-%i", $source_chr, $start - 1, $end + 1));
    next unless $iter;

    while(my $line = $iter->next) {
      chomp $line;
      my $var = $as->parse_variation($line);

      my $slice = Bio::EnsEMBL::Slice->new(
        -seq_region_name  => $var->{chr},
        -start            => $var->{start},
        -end              => $var->{end},
        -strand           => $var->{strand},
      );

      my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
        -variation_name        => $var->{variation_name},
        -seq_region_name       => $var->{chr},
        -start                 => $var->{start},
        -end                   => $var->{end},
        -slice                 => $slice,
        -strand                => $var->{strand},
        -allele_string         => $var->{allele_string},
        -is_somatic            => $var->{somatic},
        -clinical_significance => $var->{clin_sig} ? [split /,/, $var->{clin_sig}] : []
      );
      push @$variants, $vf;
    }
  } elsif (defined($as = $self->_get_AnnotationSource('Variation'))) {
    warn("RiboseqORFs plugin with non-indexed VEP cache is slow; for optimal performance, please use indexed VEP cache\n") unless $self->{slow_warning};
    $self->{slow_warning} = 1;

    # code based on AnnotationSource::Cache::Variation and AnnotationSource
    my $cache_region_size = $as->{cache_region_size};
    my ($source_chr, $min, $max, $seen, @regions) = $as->get_regions_from_coords(
      $chr, $start, $end, undef, $cache_region_size, $as->up_down_size());

    for my $region (@regions) {
      my ($c, $s) = @$region;

      my $file = $as->get_dump_file_name(
        $c,
        ($s * $cache_region_size) + 1,
        ($s + 1) * $cache_region_size
      );
      next unless -e $file;
      my $gz = gzopen($file, 'rb');

      my $line;
      while($gz->gzreadline($line)) {
        chomp $line;

        # ignore non-overlapping variants
        my $var = $as->parse_variation($line);
        next unless overlap($start, $end, $var->{start}, $var->{end});

        my $slice = Bio::EnsEMBL::Slice->new(
          -seq_region_name  => $c,
          -start            => $var->{start},
          -end              => $var->{end},
          -strand           => $var->{strand},
        );

        my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
          -variation_name        => $var->{variation_name},
          -seq_region_name       => $c,
          -start                 => $var->{start},
          -end                   => $var->{end},
          -slice                 => $slice,
          -strand                => $var->{strand},
          -allele_string         => $var->{allele_string},
          -is_somatic            => $var->{somatic},
          -clinical_significance => $var->{clin_sig} ? [split /,/, $var->{clin_sig}] : []
        );
        push @$variants, $vf;
      }
    }
  } else {
    die "ERROR: could not get variants from VEP cache";
  }
  return $variants;
}

sub _fetch_database_vars {
  # fetch variants from Ensembl API
  my ($self, $chr, $start, $end) = @_;

  my $config  = $self->{config};
  my $reg     = $config->{reg};
  my $species = $config->{species};

  $self->{sa}  ||= $reg->get_adaptor($species, 'core', 'slice');
  $self->{vfa} ||= $reg->get_adaptor($species, 'variation', 'variationfeature');

  my $slice = $self->{sa}->fetch_by_region('chromosome', $chr, $start, $end);
  next unless defined $slice;
  return $self->{vfa}->fetch_all_by_Slice($slice);
}

sub _get_paralogue_vars {
  my ($self, $chr, $start, $end, $perc_cov, $perc_pos, $all_results) = @_;

  my $vcf = $self->{vcf};
  if (defined $vcf) {
    # get variants from custom VCF file
    my @data = @{$self->get_data($chr, $start - 2, $end, $vcf)};
    for my $var (@data) {
      next unless $self->_is_clinically_significant([ $var->{clinical_significance} ]);
      my $info = $self->_prepare_vcf_info($var, $perc_cov, $perc_pos);
      $all_results = $self->_join_results($all_results, { PARALOGUE_VARIANTS => $info });
    }
  } else {
    # get Ensembl variants from mapped genomic coordinates
    my $variants;
    if ($self->{config}->{cache}) {
      $variants = $self->_fetch_cache_vars($chr, $start, $end);
    } elsif ($self->{config}->{database}) {
      $variants = $self->_fetch_database_vars($chr, $start, $end);
    } else {
      die("ERROR: cannot fetch variants from cache (no cache available?) neither from Ensembl API (database mode must be enabled)");
    }

    foreach my $var (@$variants) {
      # check clinical significance (if set)
      next unless $self->_is_clinically_significant($var->{clinical_significance});
      my $res = $self->_summarise_vf($var, $perc_cov, $perc_pos);
      $all_results = $self->_join_results($all_results, { PARALOGUE_VARIANTS => $res });
    }
  }

  if ($self->{regions}) {
    my $regions = { PARALOGUE_REGIONS => "$chr:$start:$end:$perc_cov:$perc_pos" };
    $all_results = $self->_join_results($all_results, $regions);
  }
  return $all_results;
}

## PLUGIN ----------------------------------------------------------------------

sub _get_valid_fields {
  my $selected  = shift;
  my $available = shift;

  # return all available fields when using 'all'
  return $available if $selected eq 'all';
  my @fields = split(/:/, $selected);

  # check if the selected fields exist
  my @valid;
  my @invalid;
  for my $field (@fields) {
    if ( grep { $_ eq $field } @$available ) {
      push(@valid, $field);
    } else {
      push(@invalid, $field);
    }
  }

  die "ERROR: all fields given are invalid. Available fields are:\n" .
    join(", ", @$available)."\n" unless @valid;
  warn "Paralogues plugin: WARNING: the following fields are not valid and were ignored: ",
    join(", ", @invalid), "\n" if @invalid;

  return \@valid;
}

sub new {
  my $class = shift;  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();

  my $params = $self->params_to_hash();
  my $config = $self->{config};

  # Thresholds for minimum percentage of homology similarity and coverage
  $self->{min_perc_cov} = defined $params->{min_perc_cov} ? $params->{min_perc_cov} : 0;
  $self->{min_perc_pos} = defined $params->{min_perc_pos} ? $params->{min_perc_pos} : 50;
  $self->{regions}      = defined $params->{regions}      ? $params->{regions}      : 1;

  # Prepare clinical significance parameters
  my $no_clnsig = defined $params->{clnsig} && $params->{clnsig} eq 'ignore';
  $self->{clnsig_term} = $params->{clnsig} || 'pathogenic' unless $no_clnsig;

  if (defined $self->{clnsig_term}) {
    $self->{clnsig_match} = $params->{clnsig_match} || 'partial';
    die "ERROR: clnsig_match only accepts 'exact', 'partial' or 'regex'\n"
      unless grep { $self->{clnsig_match} eq $_ } ('exact', 'partial', 'regex');
  }

  # Check information to retrieve from paralogue variants
  my $vcf = $params->{vcf};
  my @fields= ('identifier', 'alleles', 'clinical_significance'); # default
  if (defined $vcf) {
    $self->{vcf} = $vcf;
    $self->add_file($vcf);

    # get INFO fields names from VCF
    my $vcf_file = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($vcf);
    my $info = $vcf_file->get_metadata_by_pragma('INFO');
    my $info_ids = [ map { $_->{ID} } @$info ];

    @fields = @{ _get_valid_fields($params->{fields}, [@VCF_FIELDS, @$info_ids]) }
      if defined $params->{fields};

    # check if clinical significance column exists
    if (defined $self->{clnsig_term}) {
      $self->{clnsig_col} = $params->{clnsig_col};
      die "ERROR: clnsig_col must be set when using a custom VCF with clnsig\n"
        unless defined $self->{clnsig_col};

      my $filename = basename $vcf;
      die "ERROR: clnsig_col $self->{clnsig_col} not found in $filename. Available INFO fields are:\n" .
        join(", ", @$info_ids)."\n" unless grep { $self->{clnsig_col} eq $_ } @$info_ids;
    }
  } elsif (defined $params->{fields}) {
    # valid fields differ if using cache or database
    my @VAR_FIELDS = $self->{config}->{cache} ? @CACHE_FIELDS : @API_FIELDS;
    @fields = @{ _get_valid_fields($params->{fields}, \@VAR_FIELDS) };
  }
  $self->{fields} = \@fields;
  #return $self if $self->{matches};

  # Check if paralogue annotation should be downloaded
  $self->{remote} = defined $params->{paralogues} && $params->{paralogues} eq 'remote';
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
  $self->add_file($self->{paralogues});
  return $self;
}

sub feature_types {
  return ['Feature'];
}

sub get_header_info {
  my $self = shift;

  my $source;
  if (defined $self->{vcf}) {
    $source = basename $self->{vcf};
  } elsif ($self->{config}->{cache}) {
    $source = 'VEP cache';
  } elsif ($self->{config}->{database}) {
    $source = 'Ensembl API';
  }

  my $fields = join(':', @{ $self->{fields} });
  my $description = "Variants from $source in paralogue proteins (colon-separated fields: $fields)";
  my $header = { PARALOGUE_VARIANTS => $description };
  $header->{PARALOGUE_REGIONS} = 'Regions used to look up paralogue variants (colon-separated fields: chromosome:start:end:perc_cov:perc_pos)' if $self->{regions};
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
  my ($self, $vf, $perc_cov, $perc_pos) = @_;

  my @var;
  my $clnsig;
  for my $field (@{ $self->{fields} }) {
    my $info;
    if ($field eq 'identifier') {
      $info = $vf->name;
    } elsif ($field eq 'chromosome') {
      $info = $vf->seq_region_name;
    } elsif ($field eq 'start') {
      $info = $vf->seq_region_start;
    } elsif ($field eq 'end') {
      $info = $vf->seq_region_end;
    } elsif ($field eq 'strand') {
      $info = $vf->strand;
    } elsif ($field eq 'alleles') {
      $info = $vf->allele_string;
    } elsif ($field eq 'clinical_significance') {
      $info = $clnsig ||= join('/', @{$vf->get_all_clinical_significance_states});
    } elsif ($field eq 'source') {
      $info = $vf->source_name;
    } elsif ($field eq 'consequence') {
      $info = $vf->display_consequence;
    } elsif ($field eq 'perc_cov') {
      $info = $perc_cov;
    } elsif ($field eq 'perc_pos') {
      $info = $perc_pos;
    } elsif ($field eq 'gene_symbol') {
      my @symbols;
      for (@{ $vf->{slice}->get_all_Genes }) {
        push(@symbols, $_->display_xref->display_id) if $_->display_xref;
      }
      $info = join('/', @symbols);
    }
    push @var, $info || '';
  }
  return join(':', @var);
}

sub _is_clinically_significant {
  my ($self, $clnsig) = @_;

  # if no filter is set, avoid filtering
  my $filter = $self->{clnsig_term};
  return 1 unless defined $filter;

  return 0 unless any { defined } @$clnsig;

  my $res;
  my $match = $self->{clnsig_match};
  if ($match eq 'partial') {
    $filter = "\Q$filter\E";
  } elsif ($match eq 'exact') {
    $filter = "^\Q$filter\E\$";
  } elsif ($match eq 'regex') {
    # no need to do anything
  } else {
    return 0;
  }
  return grep /$filter/i, @$clnsig;
}

sub _prepare_vcf_info {
  my ($self, $data, $perc_cov, $perc_pos)= @_;
  $data->{perc_cov} = $perc_cov;
  $data->{perc_pos} = $perc_pos;

  # prepare data from selected fields
  my @res;
  for my $field (@{ $self->{fields} }) {
    my $value = defined $data->{$field} ? $data->{$field} : '';
    $value =~ s/,/ /g;
    $value =~ s/[:|]/_/g;
    push @res, $value;
  }
  return join(':', @res);
}

sub run {
  my ($self, $tva) = @_;

  # Quick checks to save time
  my $vf = $tva->variation_feature;
  my $translation_start = $tva->base_variation_feature_overlap->translation_start;
  my $translation_end   = $tva->base_variation_feature_overlap->translation_end;
  return {} unless
    defined $translation_start and defined $translation_end and
    $translation_start == $translation_end;

  return $self->_get_paralogue_vars_from_annotation($tva);
}

sub parse_data {
  my ($self, $line, $file) = @_;

  # do not parse here paralogue annotation
  return $line if defined $self->{paralogues} and $file eq $self->{paralogues};

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

  # fetch data from INFO fields
  my %data = (
    'chromosome'            => $chrom,
    'start'                 => $start,
    'identifier'            => $id,
    'alleles'               => $ref.'/'.$alt,
    'ref'                   => $ref,
    'alt'                   => $alt,
    'quality'               => $qual,
    'filter'                => $filter,
  );
  for my $field ( split /;/, $info ) {
    my ($key, $value) = split /=/, $field;
    $data{$key} = $value;
  }

  my $clnsig_col = $self->{clnsig_col};
  $data{clinical_significance} = $data{$clnsig_col} if $clnsig_col;
  return \%data;
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
