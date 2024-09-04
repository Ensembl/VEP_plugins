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

 # Same but output different fields
 ./vep -i variations.vcf --cache --plugin Paralogues,vcf=/path/to/file.vcf.gz,clnsig_col=CLNSIG,fields=identifier:alleles:CLNSIG:CLNVI:GENEINFO

 # Use a file with regions matched to paralogue variants -- fastest method;
 # download 'matches' files from https://ftp.ensembl.org/pub/current_variation/Paralogues
 ./vep -i variations.vcf --cache --plugin Paralogues,matches=Paralogues.pm_homo_sapiens_113_GRCh38_clinvar_20240107.tsv.gz,clnsig=ignore

 # Same using a 'matches' file but filtering for 'pathogenic' variants (default)
 ./vep -i variations.vcf --cache --plugin Paralogues,matches=Paralogues.pm_homo_sapiens_113_GRCh38_clinvar_20240107.tsv.gz

 # Fetch all Ensembl variants in paralogue proteins using only the Ensembl API
 # (requires database access)
 ./vep -i variations.vcf --database --plugin Paralogues,mode=remote,clnsig=ignore

=head1 DESCRIPTION

 A VEP plugin that fetches variants overlapping the genomic coordinates of amino
 acids aligned between paralogue proteins. This is useful to predict the
 pathogenicity of variants in paralogue positions.

 This plugin can determine paralogue regions for a variant based on:
   1. Pre-computed matches between genomic regions and paralogue variants.
      For this approach, either download the file calculated using ClinVar variants and respective TBI from
      https://ftp.ensembl.org/pub/current_variation/Paralogues or create such matches file yourself. Details on how
      to create such 'matches' file can be found below.
   2. Ensembl paralogue annotation. These versatile annotations can look up
      paralogue regions for all variants from any species with Ensembl
      paralogues, but take longer to process.

 After retrieving the paralogue regions, this plugin fetches variants
 overlapping those regions from one of the following sources (by this order):
   1. Custom VCF via the 'vcf' parameter
   2. VEP cache (in cache/offline mode)
   3. Ensembl API (in database mode)

 To create a 'matches' file based on a custom set of variants, run VEP using
 `--plugin Paralogues,regions=1,min_perc_cov=0,min_perc_pos=0,clnsig=ignore`
 and the `--vcf` option. Afterwards, process the output of the VEP command:
 `perl -e "use Paralogues; Paralogues::prepare_matches_file('variant_effect_output.txt')"`

 Options are passed to the plugin as key=value pairs:
   matches      : Tabix-indexed TSV file with pre-computed matches between
                  genomic regions and paralogue variants (fastest method); this
                  option is incompatible with the `paralogues` and `vcf` options

   dir          : Directory with paralogue annotation (the annotation is created
                  in this folder if the paralogue annotation files do not exist)
   paralogues   : Tabix-indexed TSV file with paralogue annotation (if the file
                  does not exist, the annotation is automatically created); if
                  set to 'remote', the annotation is fetched but not stored
   vcf          : Tabix-indexed VCF file to fetch variant information (if not
                  used, variants are fetched from VEP cache in cache/offline
                  mode or Ensembl API in database mode)

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
                    - Matches file: check column names in file header

   clnsig       : Clinical significance term to filter variants (default:
                  'pathogenic'); use 'ignore' to fetch all paralogue variants,
                  regardless of clinical significance
   clnsig_match : Type of match when filtering variants based on option
                  `clnsig`: 'partial' (default), 'exact' or 'regex'
   clnsig_col   : Column name containing clinical significance in custom VCF
                  (required with `vcf` option and if `clnsig` is not 'ignore')

   min_perc_cov : Minimum alignment percentage of the peptide associated with
                  the input variant (default: 0)
   min_perc_pos : Minimum percentage of positivity (similarity) between both
                  homologues (default: 50)
   regions      : Boolean value to return regions used to look up paralogue
                  variants (default: 1)

 The tabix utility must be installed in your path to read the paralogue
 annotation, the custom VCF file and the matches file.

=cut

package Paralogues;
@EXPORT_OK = qw(&process_data);

use strict;
use warnings;

use List::Util qw(any);
use File::Basename;
use Compress::Zlib;
use File::Spec;

use Bio::SimpleAlign;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::IO::Parser::VCF4Tabix;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

our @MATCHES_FIELDS = (
  'identifier',
  'chromosome',
  'start',
  'end',
  'alleles',
  'clinical_significance',
  'perc_cov',
  'perc_pos',
);

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

our @MATCHES_HEADER = qw/chr start end feature perc_cov perc_pos var_id var_chr var_start var_end var_ref var_alt var_feature/;

## FETCH VARIANTS --------------------------------------------------------------

=head2 _create_vf
  Arg[1]     : hash
  Description: Create variant feature from a variant hash
  Returntype : Bio::EnsEMBL::Variation::VariationFeature
  Status     : Experimental
=cut

sub _create_vf {
  my ($self, $var) = @_;

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
  return $vf;
}

=head2 _fetch_cache_vars
  Arg[1]     : $seq_region_name
  Arg[2]     : $start
  Arg[3]     : $end
  Description: Fetch variant features within a given genomic region from VEP
               cache; supports indexed (faster) and non-indexed (slower) cache
  Returntype : arrayref of Bio::EnsEMBL::Variation::VariationFeature
  Status     : Experimental
=cut

sub _fetch_cache_vars {
  my ($self, $chr, $start, $end) = @_;

  my ($as, $variants);
  if (defined($as = $self->_get_AnnotationSource('VariationTabix'))) {
    # code based on AnnotationSource::Cache::VariationTabix
    my $source_chr = $as->get_source_chr_name($chr);
    my $tabix_obj = $as->_get_tabix_obj($source_chr);
    return unless $tabix_obj;

    my $iter = $tabix_obj->query(sprintf("%s:%i-%i", $source_chr, $start - 1, $end + 1));
    return unless $iter;

    while(my $line = $iter->next) {
      chomp $line;
      my $var = $as->parse_variation($line);
      push @$variants, $self->_create_vf($var);
    }
  } elsif (defined($as = $self->_get_AnnotationSource('Variation'))) {
    warn("Using non-indexed VEP cache is slow; for optimal performance, please use indexed VEP cache\n")
      unless $self->{skip_slow_warning};
    $self->{skip_slow_warning} = 1;

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
        $var->{chr} = $c;
        next unless overlap($start, $end, $var->{start}, $var->{end});
        push @$variants, $self->_create_vf($var);
      }
    }
  } else {
    die "ERROR: could not get variants from VEP cache";
  }
  return $variants;
}

=head2 _fetch_database_vars
  Arg[1]     : $seq_region_name
  Arg[2]     : $start
  Arg[3]     : $end
  Description: Fetch variant features within a given genomic region from Ensembl
               database (requires database connection)
  Returntype : arrayref of Bio::EnsEMBL::Variation::VariationFeature
  Status     : Experimental
=cut

sub _fetch_database_vars {
  my ($self, $chr, $start, $end) = @_;

  my $config  = $self->{config};
  my $reg     = $config->{reg};
  my $species = $config->{species};

  $self->{slice_adaptor} ||= $reg->get_adaptor($species, 'core', 'slice');
  $self->{vf_adaptor}    ||= $reg->get_adaptor($species, 'variation', 'variationfeature');

  my $slice = $self->{slice_adaptor}->fetch_by_region('chromosome', $chr, $start, $end);
  next unless defined $slice;
  return $self->{vf_adaptor}->fetch_all_by_Slice($slice);
}

=head2 _fetch_database_vars
  Arg[1]     : $seq_region_name
  Arg[2]     : $start
  Arg[3]     : $end
  Description: Fetch variant features within a given genomic region from:
               - VEP cache (indexed or non-indexed) if --cache is enabled
               - Ensembl database if --database is enabled
               Throws error if --cache and --database are not enabled
  Returntype : arrayref of Bio::EnsEMBL::Variation::VariationFeature
  Status     : Experimental
=cut

sub fetch_variants {
  my ($self, $chr, $start, $end) = @_;

  my $variants;
  if ($self->config->{cache}) {
    $variants = $self->_fetch_cache_vars($chr, $start, $end);
  } elsif ($self->config->{database}) {
    $variants = $self->_fetch_database_vars($chr, $start, $end);
  } else {
    die("ERROR: cannot fetch variants from cache (no cache available?) neither from Ensembl API (database mode must be enabled)");
  }
  return $variants;
}

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
  my $transcript_id = $para_tr->stable_id || '';
  return ($chr, $start, $end, $transcript_id);
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

sub _get_paralogue_vars_from_matches {
  my ($self, $tva) = @_;
  my $vf     = $tva->variation_feature;
  my $allele = $tva->base_variation_feature->alt_alleles;

  # get transcript for this feature (skip if not available)
  return {} unless $tva->can('transcript');
  my $feat   = $tva->transcript->stable_id;

  my $file  = $self->{matches};
  my $chr   = $vf->{chr};
  my ($start, $end) = $vf->{start} < $vf->{end} ?
    ($vf->{start}, $vf->{end}) :
    ($vf->{end}, $vf->{start});

  # get paralogue variants for this region
  my @data = @{$self->get_data($chr, $start, $end, $file)};

  my $all_results = {};
  foreach (@data) {
    next unless $feat eq $_->{feature};

    my $var = $_->{var};
    $all_results = $self->_prepare_paralogue_vars_output(
      $var->{chr}, $var->{start}, $var->{end}, $var->{feature},
      $_->{perc_cov}, $_->{perc_pos}, [ $var ], $all_results);
  }
  return $all_results;
}

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

    my ($chr, $start, $end, $transcript_id) = $self->_get_paralogue_coords($tva, $aln);
    next unless defined $chr and defined $start and defined $end;

    my $variants;
    if (defined $self->{vcf}) {
      # get variants from custom VCF file
      $variants = $self->get_data($chr, $start, $end, $self->{vcf});
    } else {
      # get Ensembl variants from mapped genomic coordinates
      $variants = $self->fetch_variants($chr, $start, $end);
    }

    $all_results = $self->_prepare_paralogue_vars_output(
      $chr, $start, $end, $transcript_id, $perc_cov, $perc_pos,
      $variants, $all_results);
  }
  return $all_results;
}

sub _prepare_paralogue_vars_output {
  my ($self, $chr, $start, $end, $transcript_id, $perc_cov, $perc_pos,
      $variants, $all_results) = @_;

  my $is_rest = $self->{config}->{output_format} eq 'json' || $self->{config}->{rest};

  foreach my $var (@$variants) {
    # check clinical significance (if set)
    my $cln_sig = $var->{clinical_significance};
    $cln_sig = [ $cln_sig ] if defined $self->{vcf} or defined $self->{matches};
    next unless $self->_is_clinically_significant($cln_sig);

    my $FUN = (defined $self->{matches} or defined $self->{vcf}) ?
      '_prepare_vcf_info' : '_summarise_vf';
    my $res = $self->$FUN($var, $perc_cov, $perc_pos);

    if ($is_rest) {
      my %res_hash;
      @res_hash{ @{$self->{fields}} } = @{$res};
      $res = \%res_hash;
    } else {
      $res = join(':', @{$res});
    }
    $all_results = $self->_join_results($all_results, { PARALOGUE_VARIANTS => $res });
  }

  if ($self->{regions}) {
    my @keys   = qw/chromosome start end transcript_id perc_cov perc_pos/;
    my @values = ($chr, $start, $end, $transcript_id, $perc_cov, $perc_pos);

    my $regions;
    if ($is_rest) {
      my %regions_hash;
      @regions_hash{@keys} = @values;
      $regions = \%regions_hash;
    } else {
      $regions = join(':', @values);
    }
    $all_results = $self->_join_results($all_results, { PARALOGUE_REGIONS => $regions });
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

sub _validate_matches_file {
  my ($self, $params) = @_;
  die "ERROR: options 'matches' and 'paralogues' are incompatible\n"
    if defined $params->{paralogues};
  die "ERROR: options 'matches' and 'vcf' are incompatible\n"
    if defined $params->{vcf};
  die "ERROR: 'matches=$self->{matches}': file not found\n"
    unless -e $self->{matches};
  die "ERROR: 'matches=$self->{matches}': respective TBI file not found\n"
    unless -e $self->{matches} . ".tbi" || -e $self->{matches} . ".csi";
  return $self;
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
  if (defined $params->{matches}) {
    # File with matches between variants and their paralogues regions
    my $file = $params->{matches};
    $self->{matches} = $file;
    $self->{clnsig_col} ||= 'CLNSIG';
    $self->add_file($file);
    $self->_validate_matches_file($params);

    # Read column names from matches file
    $self->{matches_col} = `tabix -H ${file}` or die $!;
    chomp $self->{matches_col};
    $self->{matches_col} = [ split /\t/, $self->{matches_col} ];

    # Remove basic columns from user-selectable columns
    for my $elem (@MATCHES_HEADER) {
      $elem = '#' . $elem if $elem eq 'chr'; # Add hash to first column name
      $self->{matches_col} = [ grep { $elem ne $_ } @{ $self->{matches_col} } ];
    }

    my $valid_fields = [ @MATCHES_FIELDS, @{$self->{matches_col}} ];
    @fields = @{ _get_valid_fields($params->{fields}, $valid_fields) }
      if defined $params->{fields};;
  } elsif (defined $vcf) {
    $self->{vcf} = $vcf;
    $self->add_file($vcf);

    # get INFO fields names from VCF
    my $vcf_file = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($vcf);
    my $info = $vcf_file->get_metadata_by_pragma('INFO');
    my $info_ids = [ map { $_->{ID} } @$info ];

    @fields = @{ _get_valid_fields($params->{fields}, [@VCF_FIELDS, @$info_ids]) }
      if defined $params->{fields};

    # check if clinical significance column exists
    if (defined $self->{clnsig_term} || defined $params->{clnsig_col}) {
      $self->{clnsig_col} = $params->{clnsig_col};
      die "ERROR: clnsig_col must be set when using a custom VCF unless clnsig=ignore\n"
        unless defined $self->{clnsig_col} or $no_clnsig;

      my $filename = basename $vcf;
      die "ERROR: clnsig_col $self->{clnsig_col} not found in $filename. Available INFO fields are:\n" .
        join(", ", @$info_ids)."\n" unless grep { $self->{clnsig_col} eq $_ } @$info_ids;
    }

    # warn if trying to return clinical significance when the user does not input its VCF field
    warn("WARNING: clnsig_col not defined; clinical significance of paralogue variants will be empty\n")
      if !defined $self->{clnsig_col} && grep { /^clinical_significance$/ } @fields;
  } elsif (defined $params->{fields}) {
    # valid fields differ if using cache or database
    my @VAR_FIELDS = $self->{config}->{cache} ? @CACHE_FIELDS : @API_FIELDS;
    @fields = @{ _get_valid_fields($params->{fields}, \@VAR_FIELDS) };
  }
  $self->{fields} = \@fields;
  return $self if $self->{matches};

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
  if (defined $self->{matches}) {
    $source = basename $self->{matches};
  } elsif (defined $self->{vcf}) {
    $source = basename $self->{vcf};
  } elsif ($self->{config}->{cache}) {
    $source = 'VEP cache';
  } elsif ($self->{config}->{database}) {
    $source = 'Ensembl API';
  }

  my $fields = join(':', @{ $self->{fields} });
  my $description = "Variants from $source in paralogue locations (colon-separated fields: $fields)";
  my $header = { PARALOGUE_VARIANTS => $description };
  $header->{PARALOGUE_REGIONS} = 'Genomic location of paralogue regions (colon-separated fields: chromosome:start:end:transcript_id:perc_cov:perc_pos)' if $self->{regions};
  return $header;
}

sub _obj_in_arrayref {
  # Returns whether $obj is equal to any element of $arrayref
  my ($arrayref, $obj) = @_;

  for my $elem (@$arrayref) {
    if (ref $elem eq 'HASH') {
      next unless ref $obj eq 'HASH';

      # Compare if keys are the same length between the two objects
      my @keys_obj = sort keys %$obj;
      my @keys_k   = sort keys %$elem;
      next unless $#keys_k eq $#keys_obj;

      # Compare keys and values between two hashes
      my $sum = 0;
      my $total = 0;
      for (@keys_obj) {
        # sum equal values
        $sum += $obj->{$_} eq $elem->{$_} ? 1 : 0;
        $total++;
      }
      return 1 if $sum eq $total;
    } else {
      return 1 if $elem eq $obj;
    }
  }
  return 0;
}

sub _join_results {
  my $self = shift;
  my $all_results = shift;
  my $res = shift;

  # Create array of results per key
  for (keys %$res) {
    $all_results->{$_} = [] unless defined $all_results->{$_};

    if (not _obj_in_arrayref($all_results->{$_}, $res->{$_})) {
      # Avoid storing duplicate elements
      push(@{ $all_results->{$_} }, $res->{$_})
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
  return \@var;
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
  return \@res;
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

  return $self->{matches} ?
    $self->_get_paralogue_vars_from_matches($tva) :
    $self->_get_paralogue_vars_from_annotation($tva);
}

sub _parse_matches {
  my ($self, $line, $file) = @_;

  my ($chrom, $start, $end, $feature, $perc_cov, $perc_pos, $var_id, $var_chr, $var_start, $var_end, $var_ref, $var_alt, $var_feature, @INFO) = split /\t/, $line;

  my %data = (
    'chromosome' => $chrom,
    'start'      => $start,
    'end'        => $end,
    'feature'    => $feature,
    'perc_cov'   => $perc_cov,
    'perc_pos'   => $perc_pos,
    'var'        => {
      'chr'        => $var_chr,
      'chromosome' => $var_chr,
      'start'      => $var_start,
      'end'        => $var_end,
      'identifier' => $var_id,
      'alleles'    => $var_ref . "/" . $var_alt,
      'feature'    => $var_feature
    }
  );

  # add all ClinVar INFO fields to variant object
  $data{var}{ $_ } = shift @INFO for @{$self->{matches_col}};

  # Save clinical significance
  my $clnsig_col = $self->{clnsig_col};
  $data{var}{clinical_significance} = $data{var}{$clnsig_col} if $clnsig_col;

  return \%data;
}

sub _parse_vcf {
  my ($self, $line) = @_;

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
    'chromosome' => $chrom,
    'start'      => $start,
    'identifier' => $id,
    'alleles'    => $ref.'/'.$alt,
    'ref'        => $ref,
    'alt'        => $alt,
    'quality'    => $qual,
    'filter'     => $filter,
  );
  for my $field ( split /;/, $info ) {
    my ($key, $value) = split /=/, $field;
    $data{$key} = $value;
  }

  my $clnsig_col = $self->{clnsig_col};
  $data{clinical_significance} = $data{$clnsig_col} if $clnsig_col;
  return \%data;
}

sub parse_data {
  my ($self, $line, $file) = @_;

  if (defined $self->{paralogues} && $file eq $self->{paralogues}) {
    return $line;
  } elsif (defined $self->{matches} && $file eq $self->{matches}) {
    return $self->_parse_matches($line, $file);
  } else {
    return $self->_parse_vcf($line);
  }
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

sub prepare_matches_file {
  # Prepare annotation containing matches between genomic regions and paralogue
  # variant location; the file is then sorted, bgzipped and tabixed

  # Use this function on VEP output with Paralogues plugin data:
  # perl -e "use Paralogues; Paralogues::prepare_matches_file('variant_effect_output.txt')"

  my $results_file = shift;
  my $outfile = shift || dirname($results_file) . '/' . 'paralogues_matches.tsv';

  $results_file = File::Spec->rel2abs($results_file);
  die "File not found: $results_file\n" unless -e $results_file;
  my $gz = gzopen($results_file, 'r') or die "Not able to open file: $results_file\n"; 

  open(OUT, '>', $outfile) or die $!;

  my ($info_fields, $csq_cols);
  my $has_header = 0;
  my $progress = 0;
  # Split each region into new rows
  while($gz->gzreadline($_) > 0) {
    chomp;

    # read INFO columns
    if ($_ =~ /^##INFO=<ID=(.*?),/) {
      if ($1 eq 'CSQ') {
        $_ =~ /.*Format: (.*)">/;
        $csq_cols = [ split /\|/, $1 ];
      } else {
        push @$info_fields, $1;
      }
    }
    next if $_ =~ /^#/;

    die "Could not find INFO attribute for CSQ; check if VCF header is okay\n" unless $csq_cols;
    print("Processing line ", $progress, "...\n") if $progress++ % 100000 == 0;

    if (!$has_header) {
      print OUT "#", join("\t", @MATCHES_HEADER, @$info_fields), "\n";
      $has_header = 1;
    }

    my ($chr, $start, $id, $ref, $alt, $qual, $filter, $info) = split "\t", $_;
    $info = { map { split /=/ } split(/;/, $info) };
    my $csq = $info->{'CSQ'};

    for my $each_csq (split /,/, $csq) {
      my %res;
      @res{@$csq_cols} = split /\|/, $each_csq;
      if (defined $res{PARALOGUE_REGIONS} and $res{PARALOGUE_REGIONS} ne '-') {
        for my $region (split /&/, $res{PARALOGUE_REGIONS}) {
          $region =~ s/:/\t/g;
          my $end = $start + length($alt) - 1;
          my $data = join("\t", $region, $id, $chr, $start, $end, $ref, $alt, $res{Feature},
            map { $_ ne "CSQ" and defined $info->{$_} ? $info->{$_} : "NA" } @$info_fields );
          print OUT $data, "\n";
        }
      }
    }
  }
  print "Sorting file...\n";
  system "cat ${outfile} | (sed -u 1q; sort -k1,1V -k2,2n -k3,3n -k7,7V -k8,8n -k9,9n) > ${outfile}_sorted && rm ${outfile}" and die $!;

  my $outfile_gz = $outfile . '.gz';
  print "Compressing file with bgzip...\n";
  system "bgzip -c ${outfile}_sorted > $outfile_gz && rm ${outfile}_sorted"
    and die $!;

  print "Indexing file with tabix...\n";
  system "tabix -s1 -b2 -e3 $outfile_gz" and die $!;

  $gz->gzclose();
  close OUT;
}

1;
